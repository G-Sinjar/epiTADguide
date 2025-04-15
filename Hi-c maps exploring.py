# installing



# importing
import numpy as np
#import matplotlib.pyplot as plt
#import seaborn as sns
import pandas as pd
import os
import cooler


# installing the caki2 data
clr = cooler.Cooler(r"C:\Users\ghaza\Documents\ghazal\Bioinformatik_Fächer\Masterarbeit_Project\Data\caki2_experimentData\4DNFIIH3SM5N.mcool::/resolutions/25000")
print(clr.info)

# installing the ACHN data
achn_clr = cooler.Cooler(r"\Users\ghaza\Documents\ghazal\Bioinformatik_Fächer\Masterarbeit_Project\Data\ACHN_experimentData\4DNFI2DQTKBA.mcool::/resolutions/25000")



fetch_object = clr.matrix(balance=False)
# this is still not the acuall matrix. It s an object that allows fetching specific parts of the contact matrix efficiently


# get the number of contacts
import numpy as np
fetch_object = clr.matrix(balance=False, sparse=True)
total_contacts = fetch_object[:].data.sum()
print(f"Total contacts: {int(total_contacts):,}")




# fetch one region
chr = "chr13"
single_region = fetch_object.fetch(chr)
# this is a submatrix

#### explore the region ####
print(single_region) # this is a numpy array which represents values between bins
print(type(single_region))
# <class 'numpy.ndarray'>

print(single_region.shape)
# (100, 100) ->> 100 rows = 100 columns -> 100 bins --> cause 1000000/10000 = 100

print(single_region.nnz)  # Number of non-zero elements. didnt work cause it is not a sparse matrix
#df = pd.DataFrame(single_region)
#print(df)

bins = clr.bins()  # Get all bins
chr_bins = bins[bins['chrom'] == chr]  # Filter for bins in chromosome 13
print("the bins in the ",chr, "are: \n" chr_bins)



# save the file for deDoc2
import numpy as np

np.savetxt("./10kb-contact_matrix13.txt", single_region, delimiter="\t", fmt="%f")
#print("Matrix saved as 'contact_matrix13.txt' with tab delimiter.")








### Extract all interactions into a DataFrame
################# didnt work cause the cool file is huge ###########
pixels = clr.pixels()[:]  # didnt work too much for memory

# Try reading only the first few rows to check if the data can be loaded
pixels = clr.pixels(limit=1000)  # worked


# Read pixels in chunks (you can adjust chunk_size)
chunk_size = 1000000  # Adjust to a number that's appropriate for your system
pixels_iterator = clr.pixels()  # Create an iterator to go through the pixels

# Initialize an empty list to store the data
pixel_data = []

# Loop through the pixels and process them in chunks
for i, chunk in enumerate(pixels_iterator):
    if i * chunk_size >= 1000000:  # Limit reading to the first 1000000 pixels
        break
    pixel_data.append(chunk)

# Convert the data into a DataFrame
pixels_df = pd.DataFrame(pixel_data)




bins = clr.bins()[:]  # Extract bin information
bins_df = pd.DataFrame(bins)
print(bins_df.head())

# Add genomic positions instead of bin IDs
merged_df = pixels_df.merge(bins_df[['chrom', 'start']], left_on='bin1_id', right_index=True)
merged_df = merged_df.rename(columns={'chrom': 'chrom1', 'start': 'start1'})

merged_df = merged_df.merge(bins_df[['chrom', 'start']], left_on='bin2_id', right_index=True)
merged_df = merged_df.rename(columns={'chrom': 'chrom2', 'start': 'start2'})

# Keep only relevant columns
final_df = merged_df[['chrom1', 'start1', 'chrom2', 'start2', 'count']]

# Save as a TSV file
final_df.to_csv("cool_data.tsv", sep='\t', index=False)

print("TSV file saved as 'cool_data.tsv'")