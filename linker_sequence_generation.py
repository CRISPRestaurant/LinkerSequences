from Bio import SeqIO
import cooked_up_libraries.bio_string_parser as bsp
import cooked_up_libraries.dna_properties as dnap
import cooked_up_libraries.idt_writer as idtw
import json

guide_rna = "GATATCAAGAGGATTGGAAA"

direct_repeat = "gttttagagctatgctgttttgaatggtcccaaaac"
direct_repeat_hairpin_end_position = 24
direct_repeat = direct_repeat.upper()

records = list(SeqIO.parse("Plasmids/pCRCT.fasta", "fasta"))
pCRCT = records[0].seq.upper()
print(dnap.get_double_stranded_dna_molecular_weight(pCRCT))

SNR52_F = "ggctagcggtaaaggtgcgcat".upper()

GAGACC_sites = [i for i in range(len(pCRCT)) if pCRCT.startswith("GAGACC", i)]
GGTCTC_sites = [i for i in range(len(pCRCT)) if pCRCT.startswith("GGTCTC", i)]

beginning_four_nucleotides = bsp.substring(pCRCT, GAGACC_sites[0] - 2, 4, left_endpoint = False)

beginning_of_blocks = "ccaaaac".upper()
end_of_blocks = "gttttagag".upper()

blocks = beginning_of_blocks + guide_rna + direct_repeat + guide_rna + direct_repeat + end_of_blocks

guide_rna_one_forward_primer_range = [0, 45]
guide_rna_one_reverse_primer_range = [4, 49]

guide_rna_two_forward_primer_range = [46, 95]
guide_rna_two_reverse_primer_range = [50, 99]

guide_rna_three_forward_primer_range = [96, 143]
guide_rna_three_reverse_primer_range = [100, 147]

primer_info = {"FORWARD": [0, 46, 50, 48], "REVERSE": [4, 46, 50, 48]}
primers_inputted = primer_info != None

def generate_linker_sequences(guides, direct_repeat, beginning_block, end_block, minimum_length_of_linker_sequence, primer_information = "IDK, GENERATE!"):
    total_block = beginning_block

    direct_repeat_regions = []

    for guide_position in range(len(guides)):
        total_block += guides[guide_position]

        if guide_position != len(guides) - 1:
            direct_repeat_start_position = len(total_block)
            total_block += direct_repeat
            direct_repeat_end_position = len(total_block) - 1

            direct_repeat_regions.append([direct_repeat_start_position, direct_repeat_end_position])
        else:
            total_block += end_block

    forward_primer_regions = []
    reverse_primer_regions = []

    forward_primers = []
    reverse_primers = []

    if primer_information != "IDK, GENERATE!":
        forward_information = primer_information["FORWARD"]
        reverse_information = primer_information["REVERSE"]

        forward_pointer = forward_information[0]
        for forward_primer_length in range(1, len(forward_information)):
            forward_primer_position_range = [forward_pointer, forward_pointer + forward_information[forward_primer_length] - 1]
            forward_pointer = forward_primer_position_range[1] + 1

            forward_primer_regions.append(forward_primer_position_range)
            forward_primers.append(bsp.substring(total_block, forward_primer_position_range[0], forward_information[forward_primer_length], left_endpoint = True))

        reverse_pointer = reverse_information[0]
        for reverse_primer_length in range(1, len(reverse_information)):
            reverse_primer_position_range = [reverse_pointer, reverse_pointer + reverse_information[reverse_primer_length] - 1]
            reverse_pointer = reverse_primer_position_range[1] + 1

            reverse_primer_regions.append(reverse_primer_position_range)
            reverse_primers.append(dnap.getAntiSenseStrand(bsp.substring(total_block, reverse_primer_position_range[0], reverse_information[reverse_primer_length], left_endpoint = True)))
    else:
        total_space_for_primers = len(total_block) - 4
        average_length_of_each_primer = int(total_space_for_primers / len(guides))

        length_of_last_primer = total_space_for_primers - (len(guides) - 1) * average_length_of_each_primer

        forward_information = [0]
        reverse_information = [4]

        for i in range(len(guides) - 1):
            forward_information.append(average_length_of_each_primer)
            reverse_information.append(average_length_of_each_primer)

        forward_information.append(length_of_last_primer)
        reverse_information.append(length_of_last_primer)

        forward_pointer = forward_information[0]
        for forward_primer_length in range(1, len(forward_information)):
            forward_primer_position_range = [forward_pointer, forward_pointer + forward_information[forward_primer_length] - 1]
            forward_pointer = forward_primer_position_range[1] + 1

            forward_primer_regions.append(forward_primer_position_range)
            forward_primers.append(bsp.substring(total_block, forward_primer_position_range[0], forward_information[forward_primer_length], left_endpoint = True))

        reverse_pointer = reverse_information[0]
        for reverse_primer_length in range(1, len(reverse_information)):
            reverse_primer_position_range = [reverse_pointer, reverse_pointer + reverse_information[reverse_primer_length] - 1]
            reverse_pointer = reverse_primer_position_range[1] + 1

            reverse_primer_regions.append(reverse_primer_position_range)
            reverse_primers.append(dnap.getAntiSenseStrand(bsp.substring(total_block, reverse_primer_position_range[0], reverse_information[reverse_primer_length], left_endpoint = True)))

    forward_linker_sequences_before_guide = []
    reverse_linker_sequences_before_guide = []

    forward_linker_sequences_after_guide = []
    reverse_linker_sequences_after_guide = []

    # Generating previous to guide primer linker sequences for forward primers
    for i in range(1, len(forward_primer_regions)):
        beginning_position = forward_primer_regions[i][0]

        relevant_direct_repeat_region = []
        for direct_repeat_region in direct_repeat_regions:
            if beginning_position < direct_repeat_region[0]:
                break

            relevant_direct_repeat_region = direct_repeat_region

        direct_repeat_beginning_position = relevant_direct_repeat_region[0]
        length_of_direct_repeat_region_of_linker = beginning_position - direct_repeat_beginning_position

        current_linker_length = length_of_direct_repeat_region_of_linker + len(beginning_block)

        if current_linker_length >= minimum_length_of_linker_sequence + 4:
            forward_linker_sequences_before_guide.append(beginning_block + bsp.substring(total_block, direct_repeat_beginning_position, length_of_direct_repeat_region_of_linker, left_endpoint = True))
        else:
            random_addition_sequence = dnap.generateRandomCRISPRArraySequence(minimum_length_of_linker_sequence - current_linker_length)
            forward_linker_sequences_before_guide.append(beginning_block + random_addition_sequence + bsp.substring(total_block, direct_repeat_beginning_position, length_of_direct_repeat_region_of_linker, left_endpoint = True))

    # Generating previous to guide primer linker sequence for reverse primers
    for i in range(1, len(reverse_primer_regions)):
        beginning_position = reverse_primer_regions[i][0]

        relevant_direct_repeat_region = []
        for direct_repeat_region in direct_repeat_regions:
            if beginning_position < direct_repeat_region[0]:
                break

            relevant_direct_repeat_region = direct_repeat_region

        direct_repeat_beginning_position = relevant_direct_repeat_region[0]
        length_of_direct_repeat_region_of_linker = beginning_position - direct_repeat_beginning_position

        current_linker_length = length_of_direct_repeat_region_of_linker + len(beginning_block)

        if current_linker_length >= minimum_length_of_linker_sequence:
            reverse_linker_sequences_before_guide.append(dnap.getAntiSenseStrand(beginning_block[4:] + bsp.substring(total_block, direct_repeat_beginning_position, length_of_direct_repeat_region_of_linker, left_endpoint = True)))
        else:
            random_addition_sequence = dnap.generateRandomCRISPRArraySequence(minimum_length_of_linker_sequence + 4 - current_linker_length)
            reverse_linker_sequences_before_guide.append(dnap.getAntiSenseStrand(beginning_block[4:] + random_addition_sequence + bsp.substring(total_block, direct_repeat_beginning_position, length_of_direct_repeat_region_of_linker, left_endpoint = True)))

    # Generating after to guide primer linker sequence for forward primers
    for i in range(0, len(forward_primer_regions) - 1):
        end_position = forward_primer_regions[i][1]

        relevant_direct_repeat_region = []
        for direct_repeat_region in reversed(direct_repeat_regions):
            if end_position > direct_repeat_region[1]:
                break

            relevant_direct_repeat_region = direct_repeat_region

        print(end_position)
        print(relevant_direct_repeat_region)

        direct_repeat_end_position = relevant_direct_repeat_region[1]
        length_of_direct_repeat_region_of_linker = direct_repeat_end_position - end_position

        current_linker_length = length_of_direct_repeat_region_of_linker + len(end_block)

        if current_linker_length >= minimum_length_of_linker_sequence:
            forward_linker_sequences_after_guide.append(bsp.substring(total_block, direct_repeat_end_position, length_of_direct_repeat_region_of_linker, left_endpoint = False) + end_block[0:5])
        else:
            random_addition_sequence = dnap.generateRandomCRISPRArraySequence(minimum_length_of_linker_sequence + 4 - current_linker_length)
            forward_linker_sequences_after_guide.append(bsp.substring(total_block, direct_repeat_end_position, length_of_direct_repeat_region_of_linker, left_endpoint = False) + random_addition_sequence + end_block[0:5])

    # Generating after to guide primer linker sequence for reverse primers
    for i in range(0, len(reverse_primer_regions) - 1):
        end_position = reverse_primer_regions[i][1]

        relevant_direct_repeat_region = []
        for direct_repeat_region in reversed(direct_repeat_regions):
            if end_position > direct_repeat_region[1]:
                break

            relevant_direct_repeat_region = direct_repeat_region

        direct_repeat_end_position = relevant_direct_repeat_region[1]
        length_of_direct_repeat_region_of_linker = direct_repeat_end_position - end_position

        current_linker_length = length_of_direct_repeat_region_of_linker + len(end_block)

        if current_linker_length >= minimum_length_of_linker_sequence + 4:
            reverse_linker_sequences_after_guide.append(dnap.getAntiSenseStrand(bsp.substring(total_block, direct_repeat_end_position, length_of_direct_repeat_region_of_linker, left_endpoint = False) + end_block))
        else:
            random_addition_sequence = dnap.generateRandomCRISPRArraySequence(minimum_length_of_linker_sequence - current_linker_length)
            reverse_linker_sequences_after_guide.append(dnap.getAntiSenseStrand(bsp.substring(total_block, direct_repeat_end_position, length_of_direct_repeat_region_of_linker, left_endpoint = False) + random_addition_sequence + end_block))

    linker_sequences_output = {}
    for guide_number in range(len(guides)):
        if guide_number == 0:
            guide_linker_information = {
                "After": [forward_linker_sequences_after_guide[0], reverse_linker_sequences_after_guide[0]],
                "Forward": forward_primers[guide_number],
                "Reverse": reverse_primers[guide_number]
            }

            linker_sequences_output[1] = guide_linker_information
        elif guide_number == len(guides) - 1:
            guide_linker_information = {
                "Before": [forward_linker_sequences_before_guide[-1], reverse_linker_sequences_before_guide[-1]],
                "Forward": forward_primers[guide_number],
                "Reverse": reverse_primers[guide_number]
            }

            linker_sequences_output[len(guides)] = guide_linker_information
        else:
            guide_linker_information = {
                "Before": [forward_linker_sequences_before_guide[guide_number - 1], reverse_linker_sequences_before_guide[guide_number - 1]],
                "After": [forward_linker_sequences_after_guide[guide_number], reverse_linker_sequences_after_guide[guide_number]],
                "Forward": forward_primers[guide_number],
                "Reverse": reverse_primers[guide_number]
            }

            linker_sequences_output[guide_number + 1] = guide_linker_information

    return [linker_sequences_output, total_block]

number_of_guides = 3
guide_gene_names = ["ADE2"] * number_of_guides
linker_sequences = generate_linker_sequences([guide_rna] * number_of_guides, direct_repeat, beginning_of_blocks, end_of_blocks, 25, primer_information = primer_info)

total_combinations_of_sequences = []
for i in range(1, number_of_guides + 1):
    combination = [i]
    total_combinations_of_sequences.append(combination.copy())
    for j in range(i + 1, number_of_guides + 1):
        combination.append(j)
        total_combinations_of_sequences.append(combination.copy())

print(json.dumps(linker_sequences[0], indent = 4))

for combination in total_combinations_of_sequences:
    total_block = ""

    forward_strand = ""
    reverse_strand = ""

    if len(combination) == 1:
        if combination[0] == 1:
            total_block += linker_sequences[0][combination[0]]["Forward"] + linker_sequences[0][combination[0]]["After"][0]

            forward_strand += linker_sequences[0][combination[0]]["Forward"] + linker_sequences[0][combination[0]]["After"][0]
            reverse_strand += linker_sequences[0][combination[0]]["Reverse"][::-1] + linker_sequences[0][combination[0]]["After"][1][::-1]
        elif combination[0] == number_of_guides:
            total_block += linker_sequences[0][combination[0]]["Before"][0] + linker_sequences[0][combination[0]]["Forward"]

            forward_strand += linker_sequences[0][combination[0]]["Before"][0] + linker_sequences[0][combination[0]]["Forward"]
            reverse_strand += linker_sequences[0][combination[0]]["Before"][1][::-1] + linker_sequences[0][combination[0]]["Reverse"][::-1]
        else:
            total_block = linker_sequences[0][combination[0]]["Before"][0] + linker_sequences[0][combination[0]]["Forward"] + linker_sequences[0][combination[0]]["After"][0]

            forward_strand += linker_sequences[0][combination[0]]["Before"][0] + linker_sequences[0][combination[0]]["Forward"] + linker_sequences[0][combination[0]]["After"][0]
            reverse_strand += linker_sequences[0][combination[0]]["Before"][1][::-1] + linker_sequences[0][combination[0]]["Reverse"] + linker_sequences[0][combination[0]]["After"][1][::-1]
    else:
        for index in range(len(combination)):
            if index == 0:
                if combination[index] == 1:
                    total_block += linker_sequences[0][combination[index]]["Forward"]

                    forward_strand += linker_sequences[0][combination[index]]["Forward"]
                    reverse_strand += (linker_sequences[0][combination[index]]["Reverse"])[::-1]
                else:
                    total_block += linker_sequences[0][combination[index]]["Before"][0] + linker_sequences[0][combination[index]]["Forward"]

                    forward_strand += linker_sequences[0][combination[index]]["Before"][0] + linker_sequences[0][combination[index]]["Forward"]
                    reverse_strand += linker_sequences[0][combination[index]]["Before"][1][::-1] + linker_sequences[0][combination[index]]["Reverse"][::-1]
            elif index == len(combination) - 1:
                if combination[index] == number_of_guides:
                    total_block += linker_sequences[0][combination[index]]["Forward"]

                    forward_strand += linker_sequences[0][combination[index]]["Forward"]
                    reverse_strand += (linker_sequences[0][combination[index]]["Reverse"])[::-1]
                else:
                    total_block += linker_sequences[0][combination[index]]["Forward"] + linker_sequences[0][combination[index]]["After"][0]

                    forward_strand += linker_sequences[0][combination[index]]["Forward"] + linker_sequences[0][combination[index]]["After"][0]
                    reverse_strand += linker_sequences[0][combination[index]]["Reverse"][::-1] + linker_sequences[0][combination[index]]["After"][1][::-1]
            else:
                total_block += linker_sequences[0][combination[index]]["Forward"]

                forward_strand += linker_sequences[0][combination[index]]["Forward"]
                reverse_strand += (linker_sequences[0][combination[index]]["Reverse"])[::-1]
    total_block += "AGAG"

    print(combination)

    molecular_weight_of_total_construct = dnap.get_sticky_ended_dna_molecular_weight(forward_strand, reverse_strand)
    print(molecular_weight_of_total_construct)

idt_raw_input_data = []

names = []
sequences = []
scales = []
purifications = []

for i in range(1, number_of_guides + 1):
    guide_information = linker_sequences[0][i]

    if "After" in guide_information:
        names.append("LS_%sA_F" % (i))
        names.append("LS_%sA_R" % (i))

        sequences.append(guide_information["After"][0])
        sequences.append(guide_information["After"][1])

        scales.append("100nm")
        scales.append("100nm")

        purifications.append("STD")
        purifications.append("STD")

    if "Before" in guide_information:
        names.append("LS_%sB_F" % (i))
        names.append("LS_%sB_R" % (i))

        sequences.append(guide_information["Before"][0])
        sequences.append(guide_information["Before"][1])

        scales.append("100nm")
        scales.append("100nm")

        purifications.append("STD")
        purifications.append("STD")

    if not primers_inputted:
        names.append("g%s_%sF" % (guide_gene_names[i - 1], i))
        names.append("g%s_%sR" % (guide_gene_names[i - 1], i))

        sequences.append(guide_information["Forward"])
        sequences.append(guide_information["Reverse"])

        scales.append("100nm")
        scales.append("100nm")

        purifications.append("STD")
        purifications.append("STD")

idt_raw_input_data = [names, sequences, scales, purifications]
idtw.write_sequences_to_idt_form(idt_raw_input_data, "IDT Forms/Linker Sequences 03-30-2021.xlsx")
