#
import pandas as pd
import os

""" Set paths """
current_path = os.getcwd()
script_path = os.path.join(current_path, "scripts")

col_type = {"structure":str,
            "chain":str,
            "num_aa":str,
            "substrate":str}

def split_SS(STR_values):
    unstructured = ["U", "S", "T", "B"]
    structured = ["E", "H", "G", "I"]
    parts_STR = []
    local_part = ''
    for index, value in enumerate(STR_values):
        if index != len(STR_values)-1:
            if (value in unstructured and STR_values[index+1] in unstructured) or (value in structured and STR_values[index+1] in structured):
                local_part += value
            else:
                local_part += value
                parts_STR.append(local_part)
                local_part = ''
        else:
            if (value in unstructured and STR_values[index-1] in unstructured) or (value in structured and STR_values[index-1] in structured):
                local_part += value
                parts_STR.append(local_part)
            else:
                 parts_STR.append(value)
    return parts_STR

def process_edge(parts_STR):
    bends = ["SS", "TT", "ST", "TS", "US", "UT", "SU", "TU", "BS", "SB", "BT", "TB"]
    edge = ''

    # part 1
    part_1 = parts_STR[0]
#    print(f"Part 1: {part_1}")
    bend_counts = 0
    if len(part_1) % 2 == 0:
        for index in range(0, len(part_1)):
            if index % 2 == 0:
                edge += part_1[index] + part_1[index+1]
            #    print(part_1[index] + part_1[index+1], end=' ')
                if part_1[index] + part_1[index+1] in bends:
                    bend_counts += 1

            if bend_counts == 3:
            #    print(f"\nBend counts: {bend_counts}")
            #    print(f"Edge: {edge}")
                return edge
        #input("#2")
    else:
        if len(part_1) == 1:
            edge = part_1
        else:
            for index in range(0, len(part_1)-1):
                if index % 2 == 0:
                    edge += part_1[index] + part_1[index+1]
                #    print(part_1[index] + part_1[index+1], end=' ')
                    if part_1[index] + part_1[index+1] in bends:
                        bend_counts += 1
                elif index == len(part_1)-2:
                    edge += part_1[index+1]
                #    print(part_1[index] + part_1[index+1], end=' ')
                    if part_1[index] + part_1[index+1] in bends:
                        bend_counts += 1

                if bend_counts == 3:
                #    print(f"\nBend counts: {bend_counts}")
                #    print(f"Edge: {edge}")
                    return edge
        #input("#1")

    #part 2
    #print(f"Bend counts: {bend_counts}")

    part_2 = parts_STR[1]
    part_3 = parts_STR[2]

    if bend_counts < 3:
        #print(f"Part 2: {part_2}")
        #print(f"Part 3: {part_3}")
        if len(set(part_2)) > 1:
            edge += part_2
        else:
            # part 3
            if len(part_3) == 1:
                #edge += part_2 + part_3 + part_4 ### ??? Split part 4 or not
                part_4 = parts_STR[3]
                #print(f"Part 4: {part_4}")
                if len(set(part_4)) > 1:
                    if ("H" in set(part_4) and "I" in set(part_4)) or ("H" in set(part_4) and "G" in set(part_4)) or ("G" in set(part_4) and "I" in set(part_4)):
                        edge += part_2 + part_3 + part_4
                    else:
                        for i, l in enumerate(part_4):
                            if l != part_4[i+1]:
                                edge += part_2 + part_3 + part_4[:i+1]
                                break
                else:
                    edge += part_2 + part_3 + part_4
            elif len(part_3) == 2:
                if part_3 in bends:
                    edge += part_2 + part_3
                else:
                    part_4 = parts_STR[3]
                    #print(f"Part 4: {part_4}")
                    if len(set(part_4)) > 1:
                        if ("H" in set(part_4) and "I" in set(part_4)) or ("H" in set(part_4) and "G" in set(part_4)) or ("G" in set(part_4) and "I" in set(part_4)):
                            edge += part_2 + part_3 + part_4
                        else:
                            for i, l in enumerate(part_4):
                                if l != part_4[i+1]:
                                    edge += part_2 + part_3 + part_4[:i+1]
                                    break
                    else:
                        edge += part_2 + part_3 + part_4
            else:
                if part_3[:2] in bends:
                    edge += part_2 + part_3[:2]
                elif part_3[:3] in ["UUU", "UUS", "UUT", "UUB", "BUS", "BUT", "BUU", "UBS", "UBT", "UBU"]:
                    edge += part_2 + part_3[:3]
                else:
                    input("*")
    #print(f"Edge: {edge}")
    #input()
    return edge

def main():
    feature_files = [i for i in os.listdir(script_path) if '.csv' in i]

    for feature_file in feature_files:
        feature_df = pd.read_csv(os.path.join(script_path, feature_file), dtype=col_type)

        parts_STR = split_SS(feature_df["init_SS_type"].values)
        #print(parts_STR)
        if len(parts_STR) <= 3:
            input("!")
        N_terminal_edge = process_edge(parts_STR)
        C_terminal_edge = process_edge([i[::-1] for i in parts_STR[::-1]])

        edge_part_feature = [1]*len(N_terminal_edge) + [0]*(len(feature_df) - len(N_terminal_edge) - len(C_terminal_edge)) + [1]*len(C_terminal_edge)
        feature_df["is_edge"] = edge_part_feature

        feature_df.to_csv(os.path.join(script_path, feature_file), index=False)

if __name__ == "__main__":
    main()
