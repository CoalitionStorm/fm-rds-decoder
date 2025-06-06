import numpy as np
from scipy.signal import resample_poly, firwin
from rtlsdr import RtlSdr

# Constants
syndrome = [383, 14, 303, 663, 748]
offset_pos = [0, 1, 2, 3, 2]
offset_word = [252, 408, 360, 436, 848]
pty_table = [["Undefined", "Undefined"], ["News", "News"], ["Current Affairs", "Information"],
             ["Information", "Sports"], ["Sport", "Talk"], ["Education", "Rock"], ["Drama", "Classic Rock"],
             ["Culture", "Adult Hits"], ["Science", "Soft Rock"], ["Varied", "Top 40"], ["Pop Music", "Country"],
             ["Rock Music", "Oldies"], ["Easy Listening", "Soft"], ["Light Classical", "Nostalgia"],
             ["Serious Classical", "Jazz"], ["Other Music", "Classical"], ["Weather", "Rhythm & Blues"],
             ["Finance", "Soft Rhythm & Blues"], ["Children’s Programmes", "Language"],
             ["Social Affairs", "Religious Music"], ["Religion", "Religious Talk"], ["Phone-In", "Personality"],
             ["Travel", "Public"], ["Leisure", "College"], ["Jazz Music", "Spanish Talk"],
             ["Country Music", "Spanish Music"], ["National Music", "Hip Hop"], ["Oldies Music", "Unassigned"],
             ["Folk Music", "Unassigned"], ["Documentary", "Weather"], ["Alarm Test", "Emergency Test"],
             ["Alarm", "Emergency"]]
coverage_area_codes = ["Local", "International", "National", "Supra-regional",
                       "Regional 1", "Regional 2", "Regional 3", "Regional 4", "Regional 5",
                       "Regional 6", "Regional 7", "Regional 8", "Regional 9", "Regional 10",
                       "Regional 11", "Regional 12"]
pty_locale = 1

# Globals
radiotext_AB_flag = 0
radiotext = [' '] * 65
first_time = True

last_text = ""

def parse_rds(group):
    global radiotext_AB_flag, radiotext, first_time, last_text

    group_0 = group[1] | (group[0] << 8)
    group_1 = group[3] | (group[2] << 8)
    group_2 = group[5] | (group[4] << 8)
    group_3 = group[7] | (group[6] << 8)

    group_type = (group_1 >> 12) & 0xf
    AB = (group_1 >> 11) & 0x1
    program_type = (group_1 >> 5) & 0x1f
    pty = pty_table[program_type][pty_locale]

    if first_time:
        pi_area_coverage = (group_0 >> 8) & 0xf
        coverage_area = coverage_area_codes[pi_area_coverage]
        pi_program_reference_number = group_0 & 0xff
        print("PTY:", pty)
        print("program:", pi_program_reference_number)
        print("coverage_area:", coverage_area)
        first_time = False

    if group_type == 2:
        if radiotext_AB_flag != ((group_1 >> 4) & 0x01):
            radiotext = [' '] * 65
        radiotext_AB_flag = (group_1 >> 4) & 0x01
        text_segment_address_code = group_1 & 0x0f
        if AB:
            radiotext[text_segment_address_code * 2] = chr((group_3 >> 8) & 0xff)
            radiotext[text_segment_address_code * 2 + 1] = chr(group_3 & 0xff)
        else:
            radiotext[text_segment_address_code * 4] = chr((group_2 >> 8) & 0xff)
            radiotext[text_segment_address_code * 4 + 1] = chr(group_2 & 0xff)
            radiotext[text_segment_address_code * 4 + 2] = chr((group_3 >> 8) & 0xff)
            radiotext[text_segment_address_code * 4 + 3] = chr(group_3 & 0xff)

        current_text = ''.join(radiotext).strip()
        if len(current_text) > 5 and current_text != last_text:
            print("Radiotext:", current_text)
            last_text = current_text
 

def decode_rds(bits):
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

    synced = False
    presync = False
    wrong_blocks_counter = 0
    blocks_counter = 0
    group_good_blocks_counter = 0
    reg = 0
    lastseen_offset_counter = 0
    lastseen_offset = 0
    block_number = 0
    block_bit_counter = 0
    group_assembly_started = False
    bytes_out = []

    for i in range(len(bits)):
        reg = int(((reg << 1) | int(bits[i])) & 0x3FFFFFF)
        if not synced:
            reg_syndrome = calc_syndrome(reg, 26)
            for j in range(5):
                if reg_syndrome == syndrome[j]:
                    if not presync:
                        lastseen_offset = j
                        lastseen_offset_counter = i
                        presync = True
                    else:
                        block_distance = (offset_pos[j] - offset_pos[lastseen_offset]) % 4
                        if block_distance * 26 == (i - lastseen_offset_counter):
                            synced = True
                            block_number = (j + 1) % 4
                            group_assembly_started = False
                    break
        else:
            if block_bit_counter < 25:
                block_bit_counter += 1
            else:
                good_block = False
                dataword = (reg >> 10) & 0xffff
                checkword = reg & 0x3ff
                block_received_crc = checkword ^ offset_word[block_number]
                block_calculated_crc = calc_syndrome(dataword, 16)
                if block_received_crc == block_calculated_crc:
                    good_block = True
                elif block_number == 2:  # special check for block C
                    block_received_crc = checkword ^ offset_word[4]
                    if block_received_crc == block_calculated_crc:
                        good_block = True
                else:
                    wrong_blocks_counter += 1

                if block_number == 0 and good_block:
                    group_assembly_started = True
                    group_good_blocks_counter = 1
                    group = bytearray(8)
                if group_assembly_started:
                    if good_block:
                        group[block_number*2] = int((dataword >> 8) & 255)
                        group[block_number*2+1] = int(dataword & 255)
                        group_good_blocks_counter += 1
                        if group_good_blocks_counter == 5:
                            bytes_out.append(group)
                    else:
                        group_assembly_started = False

                block_bit_counter = 0
                block_number = (block_number + 1) % 4
                blocks_counter += 1
                if blocks_counter == 50:
                    if wrong_blocks_counter > 35:
                        synced = False
                        presync = False
                    blocks_counter = 0
                    wrong_blocks_counter = 0
    return bytes_out

def symbol_sync(samples, sps=16, mu_init=0.01):
    samples_interpolated = resample_poly(samples, 32, 1)
    mu = mu_init
    out = np.zeros(len(samples), dtype=np.complex64)
    i_in, i_out = 0, 2
    while i_out < len(samples) and i_in+32 < len(samples):
        out[i_out] = samples_interpolated[i_in*32 + int(mu*32)]
        x = (np.real(out[i_out]) > 0) + 1j * (np.imag(out[i_out]) > 0)
        y = out[i_out]
        mm_val = np.real((y - out[i_out-2]) * np.conj(x - out[i_out-2]))
        mu += sps + 0.01 * mm_val
        i_in += int(np.floor(mu))
        mu -= np.floor(mu)
        i_out += 1
    return out[2:i_out]

def costas_loop(samples, alpha=8.0, beta=0.02):
    N = len(samples)
    phase = 0
    freq = 0
    out = np.zeros(N, dtype=np.complex64)
    for i in range(N):
        out[i] = samples[i] * np.exp(-1j * phase)
        error = np.real(out[i]) * np.imag(out[i])
        freq += beta * error
        phase += freq + alpha * error
        phase %= 2 * np.pi
    return out

# SDR setup
sdr = RtlSdr()
sdr.sample_rate = 250e3
sdr.center_freq = 104.1e6  # Change to a known RDS station
sdr.gain = 'auto'

print("Starting FM RDS decoding...")

while True:
    samples = sdr.read_samples(256*1024).astype(np.complex64)

    x = 0.5 * np.angle(samples[0:-1] * np.conj(samples[1:])).astype(np.complex64)
    t = np.arange(len(x)) / sdr.sample_rate
    x *= np.exp(2j * np.pi * -57e3 * t).astype(np.complex64)

    # Bandpass filter around 57kHz
    x = np.convolve(x, firwin(101, 7.5e3, fs=sdr.sample_rate), 'valid')

    # Decimate and resample
    x = x[::10]
    x = resample_poly(x, 19, 25)

    x = symbol_sync(x)
    x = costas_loop(x)

    bits = (np.real(x) > 0).astype(np.uint8)
    bits = (bits[1:] - bits[:-1]) % 2  # Differential decoding

    groups = decode_rds(bits)
    for group in groups:
        parse_rds(group)

