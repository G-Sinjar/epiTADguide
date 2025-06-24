# importing libraries
import os
import numpy as np
import cooler
import re


def process_contact_matrix(tissue, chrom, resolution_bp, path_mcool_file):
    try:
        # ✅ Normalize chromosome name
        chrom = chrom.lower()

        # ✅ Check resolution_bp
        resolution_str = str(resolution_bp).lower().strip()
        # Regex for accepted formats:
        # - digits only, e.g. "25000"
        # - digits followed by 'k', e.g. "25k"
        pattern = r'^\d+k?$'  # matches digits optionally followed by one 'k'

        if not re.match(pattern, resolution_str):
            return {
                "success": False,
                "message": f"❌ Invalid resolution format: '{resolution_bp}'\nAllowed formats are either digits only (e.g. '25000') or digits followed by 'k' (e.g. '25k').",
                "error": "Invalid resolution input"
            }

        # Now parse the input
        if resolution_str.endswith('k'):
            resolution_bp = int(resolution_str[:-1]) * 1000
        else:
            resolution_bp = int(resolution_str)

        # Derive bin size from cleaned resolution
        binsize = str(resolution_bp // 1000)
        results_folder = f"{os.path.dirname(path_mcool_file)}/TADs_{tissue}"
        os.makedirs(f"{results_folder}/contact_matrix", exist_ok=True)


        # ✅ Check if .mcool file exists
        mcool_full_path = f"./{path_mcool_file}.mcool"
        if not os.path.exists(mcool_full_path):
            return {
                "success": False,
                "message": f"❌ File not found: {mcool_full_path}\nPlease check the path and try again.",
                "error": "Missing .mcool file"
            }

        # ✅ Load Cooler object
        clr_path = f"{mcool_full_path}::/resolutions/{resolution_bp}"
        clr = cooler.Cooler(clr_path)

        # ✅ Fetch matrix for chromosome
        fetch_object = clr.matrix(balance=False)
        submatrix = fetch_object.fetch(chrom)

        # ✅ Save to output file
        output_path = f"{results_folder}/contact_matrix/{tissue}_{chrom}_{binsize}kb_contact_matrix.txt"
        np.savetxt(output_path, submatrix, delimiter="\t", fmt="%f")

        # ✅ Get total contact count
        sparse_matrix = clr.matrix(balance=False, sparse=True)
        total_contacts = sparse_matrix[:].data.sum()

        return {
            "success": True,
            "output_file": output_path,
            "message": f"✔️ Contact matrix saved at: {output_path}",
            "total_contacts": int(total_contacts)
        }

    except Exception as e:
        return {
            "success": False,
            "message": "❌ Failed to process contact matrix due to an unexpected error.",
            "error": str(e)
        }


# ## trying the function
# result = process_contact_matrix(
#     tissue="CAKI2",
#     chrom="Chr5",
#     resolution_bp = "25k",
#     path_mcool_file="./test/4DNFIIH3SM5N"
# )
#
# if result["success"]:
#     print(result["message"])
#     print(f"Total contacts: {result['total_contacts']:,}")
# else:
#     print(result["message"])
#     print("Error:", result["error"])









######## this part is now written in r and integrated in the tadcalling shiny app
# ### deDock2 integrate
# import subprocess
#
# # Your dynamic variables
# jar_path = "./deDoc2-main/deDoc2.jar"
# input_matrix = f"./{caki2_results_folder}/contact_matrix/{tissue}_{chr}_{binsize}kb_contact_matrix.txt"
# output_tads = f"./{caki2_results_folder}/2_TADcaller_results/{tissue}_{chr}_{binsize}kb_TADlines"
#
# # Command as list of arguments (preferred)
# command = ["java", "-Xms16g","-jar", jar_path, "-inputfile", input_matrix,"-binsize", binsize,"-outputfile", output_tads]
# print("This precess can take a few minutes...")
# # Run the command
# deDock2_result = subprocess.run(command, capture_output=True, text=True)
#
# # Optional: check or print result
# if deDock2_result.returncode == 0:
#     print("TADcaller deDock2 ran successfully.")
# else:
#     print("Error running TADcaller deDock2:", deDock2_result.stderr)
# ###











