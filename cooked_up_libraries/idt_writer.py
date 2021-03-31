import numpy as np
import pandas as pd

def write_sequences_to_idt_form(sequences, filepath):
    idt_dataframe = pd.DataFrame({"Name": sequences[0], "Sequence": sequences[1], "Scale": sequences[2], "Purification": sequences[3]})

    idt_dataframe.to_excel(filepath, index = False)
