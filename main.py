import numpy as np
from scipy.signal import resample_poly, firwin, bilinear, lfilter
import matplotlib.pyplot as plt

#Read in signal
x = np.fromfile('C:\\Users\\Soham\\Downloads\\fm_rds_250k_1Msamples.iq', dtype=np.complex64)
sample_rate = 250e3
center_freq = 99.5e6

from scipy.io import wavfile

#Demodulation
x = np.diff(np.unwrap(np.angle(x)))

#De-emphasis filter
bz, az = bilinear(1, [75e-6, 1], fs=sample_rate)
x = lfilter(bz, az, x)

#decimate by 6 to get mono audio
x = x[::6]
sample_rate_audio = sample_rate/6

#normalize volume so its between -1 and +1
x /= np.max(np.abs(x))

#some machines want int16s
x *= 32767
x = x.astype(np.int16)

#Save to wav file, to open in Audacity
wavfile.write('fm.wav', int(sample_rate_audio), x)



#Quadrature Demodulation
x = 0.5 * np.angle(x[0:-1] * np.conj(x[1:])) 

#Frequency Shift
N = len(x)
f_o = -57e3 # amount to shift by
t = np.arange(N)/sample_rate # time vector
x = x * np.exp(2j*np.pi*f_o*t)

#Filter to Isolate RDS
taps = firwin(numtaps = 101, cutoff = 7.5e3, fs = sample_rate)
x = np.convolve(x, taps, 'valid')

# Preventing Aliasing through Decimation
x = x[::10]
sample_rate = 25e3

#Resample to 19 kHz
x = resample_poly(x, 19, 25)
sample_rate = 19e3

#Symbol_sync
samples = x
samples_interpolated = resample_poly(samples, 32 , 1) #32 is interpolation factor
sps = 16
mu = 0.01
out = np.zeros(len(samples) + 10, dtype = np.complex64)
out_rail = np.zeros(len(samples) + 10, dtype = np.complex64)
i_in = 0
i_out = 2
while i_out < len(samples) and i_in*32 < len(samples):
    out[i_out] = samples_interpolated[i_in*32 + int(mu*32)]
    out_rail[i_out] = int(np.real(i_out) > 0) + 1j*int(np.imag(out[i_out-1]))
    x = (out_rail[i_out] - out_rail[i_out-2]) * np.conj(out[i_out-1])
    y = (out[i_out] - out[i_out-2]) * np.conj(out_rail[i_out-1])
    mm_val = np.real(y - x)
    mu += sps + 0.01*mm_val
    i_in += int(np.floor(mu))
    mu = mu - np.floor(mu)
    i_out += 1
x = out[2:i_out]

#Fine Frequency Sync
samples = x
N = len(samples)
phase = 0
freq = 0
#These two parameters make feedback loop faster or slower
alpha = 8.0
beta = 0.2
out = np.zeros(N, dtype=np.complex64)
freq_log = []
for i in range(N):
    out[i] = samples[i] * np.exp(-1j*phase) #adjust input sample by inverse of phase offset
    error = np.real(out[i]) * np.imag(out[i])

    #Used to advance loop
    freq += (beta*error)
    freq_log.append(freq * sample_rate / (2*np.pi))
    phase += freq + (alpha * error)

    # Phase wraps around 0 or 2pi
    while phase >= 2*np.pi:
        phase -= 2*np.pi
    while phase < 0:
        phase += 2*np.pi
    
x = out

#Demodulating BPSK
bits = (np.real(x) > 0).astype(int)

#Differential Decoding
bits = (bits[1:] - bits[0:-1]) % 2
bits = bits.astype(np.uint8)

#RDS Decoding

#Constants
syndrome = [383, 14, 303, 663, 748]
offset_pos = [0, 1, 2, 3, 2]
offset_word = [252, 408, 360, 436, 848]

def calc_syndrome(x, mlen):
    reg = 0
    plen = 10
    for ii in range(mlen, 0, -1):
        reg = (reg << 1) | ((x >> (ii-1)) & 0x01)
        if (reg & (1 << plen)):
            reg = reg ^ 0x5B9
    for ii in range(plen, 0, -1):
        reg = reg << 1
        if (reg & (1 << plen)):
            reg = reg ^ 0x5B9
    return reg & ((1 << plen) - 1)

#Initalize all the working vars we will need during loop
synced = False
presync = False

wrong_blocks_counter = 0
blocks_counter = 0
group_good_blocks_counter = 0

reg = np.uint32(0)
lastseen_offset_counter = 0
lastseen_offset = 0

bytes_out = []
for i in range(len(bits)):
    reg = np.bitwise_or(np.left_shift(reg, 1), bits[i])
    if not synced:
        reg_syndrome = calc_syndrome(reg, 26)
        for j in range(5):
            if reg_syndrome == syndrome[j]:
                if not presync:
                    lastseen_offset = j
                    lastseen_offset_counter = i
                    presync = True
                else:
                    if offset_pos[lastseen_offset] >= offset_pos[j]:
                        block_distance = offset_pos[j] + 4 - offset_pos[lastseen_offset]
                    else:
                        block_distance = offset_pos[j] - offset_pos[lastseen_offset]
                    if (block_distance*26) != (i - lastseen_offset_counter):
                        presync = False
                    else:
                        print('Sync Stated Detected')
                        wrong_blocks_counter = 0
                        blocks_counter = 0
                        block_bit_counter = 0
                        block_number = (j+1) % 4
                        group_assembly_started = False
                        synced = True
            break

        else:
            if block_bit_counter < 25:
                block_bit_counter += 1
            else:
                good_block = False
                dataword   = (reg >> 10) & 0xffff
                block_calculated_crc = calc_syndrome(dataword, 16)
                checkword = reg & 0x3ff
                if block_number == 2:
                    block_received_crc = checkword ^ offset_word[block_number]
                    if (block_received_crc == block_calculated_crc):
                        good_block = True
                    else:
                        wrong_blocks_counter += 1
                        good_block = False

        if block_number == 0 and good_block:
            group_assembly_started = True
            group_good_blocks_counter = 1
            group = bytearray(8)
        if group_assembly_started:
            if not good_block:
                group_assembly_started = False
            else:
                group[block_number*2] = (dataword >> 8) & 255
                group[block_number*2+1] = dataword & 255
                group_good_blocks_counter += 1

            if group_good_blocks_counter == 5:
                bytes_out.append(group)
        
        block_bit_counter = 0
        block_number = (block_number + 1) % 4
        blocks_counter += 1
        if blocks_counter == 50:
            if wrong_blocks_counter > 35:
                print("Lost Sync (Got ", wrong_blocks_counter, " bad blocks on ", blocks_counter, " total)")
                synced = False
                presync = False
            else: 
                print("Still Sync-ed (Got ", wrong_blocks_counter, " bad blocks on ", blocks_counter, "total )")
            blocks_counter = 0
            wrong_blocks_counter = 0

pty_table = [["Undefined",             "Undefined"],
             ["News",                  "News"],
             ["Current Affairs",       "Information"],
             ["Information",           "Sports"],
             ["Sport",                 "Talk"],
             ["Education",             "Rock"],
             ["Drama",                 "Classic Rock"],
             ["Culture",               "Adult Hits"],
             ["Science",               "Soft Rock"],
             ["Varied",                "Top 40"],
             ["Pop Music",             "Country"],
             ["Rock Music",            "Oldies"],
             ["Easy Listening",        "Soft"],
             ["Light Classical",       "Nostalgia"],
             ["Serious Classical",     "Jazz"],
             ["Other Music",           "Classical"],
             ["Weather",               "Rhythm & Blues"],
             ["Finance",               "Soft Rhythm & Blues"],
             ["Children’s Programmes", "Language"],
             ["Social Affairs",        "Religious Music"],
             ["Religion",              "Religious Talk"],
             ["Phone-In",              "Personality"],
             ["Travel",                "Public"],
             ["Leisure",               "College"],
             ["Jazz Music",            "Spanish Talk"],
             ["Country Music",         "Spanish Music"],
             ["National Music",        "Hip Hop"],
             ["Oldies Music",          "Unassigned"],
             ["Folk Music",            "Unassigned"],
             ["Documentary",           "Weather"],
             ["Alarm Test",            "Emergency Test"],
             ["Alarm",                 "Emergency"]]

pty_locale = 1

coverage_area_codes = ["Local",
                       "International",
                       "National",
                       "Supra-regional",
                       "Regional 1",
                       "Regional 2",
                       "Regional 3",
                       "Regional 4",
                       "Regional 5",
                       "Regional 6",
                       "Regional 7",
                       "Regional 8",
                       "Regional 9",
                       "Regional 10",
                       "Regional 11",
                       "Regional 12"]

radiotext_AB_flag = 0
radiotext = [' ']*65
first_time = True

for group in bytes_out:
    group_0 = group[1] | (group[0] << 8)
    group_1 = group[3] | (group[2] << 8)
    group_2 = group[5] | (group[4] << 8)
    group_3 = group[7] | (group[6] << 8)

    group_type = (group_1 >> 12)
    AB = (group_1 >> 11)

    program_identification = group_0

    program_type = (group_1 >> 5) & 0x1f
    pty = pty_table[program_type][pty_locale]

    pi_area_coverage = (program_identification >> 8) & 0xf
    coverage_area = coverage_area_codes[pi_area_coverage]

    pi_program_reference_number = program_identification & 0xff

    if first_time:
        print("PTY:", pty)
        print("program:", pi_program_reference_number)
        print("coverage_area:", coverage_area)
        first_time = False

    if group_type == 2:
        if radiotext_AB_flag != ((group_1 >> 4) & 0x01):
            radiotext = [' ']*65
        radiotext_AB_flag = (group_1 >> 4) & 0x01
        text_segment_address_code = group_1 & 0x0f

        if AB:
            radiotext[text_segment_address_code * 2    ] = chr((group_3 >> 8) & 0xff)
            radiotext[text_segment_address_code * 2 + 1] = chr(group_3        & 0xff)
        else:
            radiotext[text_segment_address_code *4     ] = chr((group_2 >> 8) & 0xff)
            radiotext[text_segment_address_code * 4 + 1] = chr(group_2        & 0xff)
            radiotext[text_segment_address_code * 4 + 2] = chr((group_3 >> 8) & 0xff)
            radiotext[text_segment_address_code * 4 + 3] = chr(group_3        & 0xff)
        print(''.join(radiotext))
    else:
        pass


