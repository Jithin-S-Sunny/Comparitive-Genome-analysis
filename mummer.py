import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
alignment_file = "alignment_summary.txt"  
ref_start, ref_end, strain_start, strain_end, identity = [], [], [], [], []
with open(alignment_file, "r") as file:
    lines = file.readlines()
for line in lines[6:]:  
    fields = line.strip().split('|')
    if len(fields) >= 4:
        try:
            ref_s, ref_e = map(int, fields[0].split())
            strain_s, strain_e = map(int, fields[1].split())
            idy = float(fields[3].strip())

            ref_start.append(ref_s)
            ref_end.append(ref_e)
            strain_start.append(strain_s)
            strain_end.append(strain_e)
            identity.append(idy)
        except ValueError:
            continue  
df_plot = pd.DataFrame({
    "Ref_Start": ref_start,
    "Ref_End": ref_end,
    "Strain_Start": strain_start,
    "Strain_End": strain_end,
    "Identity": identity
})
plt.figure(figsize=(10, 6))

sns.scatterplot(x=df_plot["Ref_Start"], y=df_plot["Strain_Start"], s=50, color="blue", label="Start Positions")
sns.scatterplot(x=df_plot["Ref_End"], y=df_plot["Strain_End"], s=50, color="red", label="End Positions")
for i, row in df_plot.iterrows():
    plt.plot([row["Ref_Start"], row["Ref_End"]], [row["Strain_Start"], row["Strain_End"]],
             linestyle="-", color="black", alpha=0.6)
plt.xlabel("Reference Genome Position (AE004091.2)")
plt.ylabel("Strain Genome Position (lcl|contig_1)")
plt.title("MUMmer-Style Dot Plot of Structural Variants")
plt.legend()
plt.grid(True)
plt.show()

