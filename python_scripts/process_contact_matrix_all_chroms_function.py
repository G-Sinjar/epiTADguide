import os
import numpy as np
import cooler


#------------------------------------------------
def process_contact_matrix_all_chroms(tissue, resolution_kb, path_mcool_file):
    """
    Args:
        tissue (str): The tissue name (e.g., "CAKI2").
        resolution_kb (int): The resolution in kilobases (e.g., 25).
        path_mcool_file (str): The path to the .mcool file. Can be:
                                 - Full path with .mcool extension (e.g., "/data/file.mcool")
                                 - Path without .mcool extension (e.g., "/data/file")
    Returns:
        dict: A dictionary containing:
              - "success" (bool): True if all chromosomes were processed successfully, False otherwise.
              - "chromosomes_processed" (list): A list of dictionaries, each with:
                    - "chrom" (str): Chromosome name.
                    - "output_file" (str, optional): Path to the saved contact matrix file for this chromosome.
                    - "message" (str): A status or error message for this chromosome.
                    - "error" (str, optional): A brief error description if not successful for this chromosome.
              - "total_contacts_overall_mcool" (int, optional): Sum of all contacts in the *entire* .mcool file.
              - "overall_message" (str): A summary message for the entire operation.
    """
    overall_success = True
    processed_chromosomes_info = []
    total_contacts_overall = 0

    try:
        # --- Input Normalization and Validation ---
        try:
            resolution_kb = int(resolution_kb)
            resolution_bp = resolution_kb * 1000
        except ValueError:
            return {
                "success": False,
                "overall_message": f"❌ Invalid resolution input: '{resolution_kb}'. Resolution must be an integer in kilobases (e.g., 25).",
                "chromosomes_processed": [],
                "total_contacts_overall_mcool": 0
            }

        mcool_full_path = path_mcool_file
        if not mcool_full_path.endswith(".mcool"):
            mcool_full_path = f"{path_mcool_file}.mcool"

        mcool_full_path_abs = os.path.abspath(mcool_full_path)

        if not os.path.exists(mcool_full_path_abs):
            return {
                "success": False,
                "overall_message": f"❌ File not found: '{mcool_full_path_abs}'. Please check the path and filename (with or without .mcool extension) and try again.",
                "chromosomes_processed": [],
                "total_contacts_overall_mcool": 0
            }

        if os.path.isdir(mcool_full_path_abs):
            return {
                "success": False,
                "overall_message": f"❌ Input path resolves to a directory: '{mcool_full_path_abs}'. \nPlease provide the full path to the .mcool file, including its name (e.g., /data/my_file.mcool).",
                "chromosomes_processed": [],
                "total_contacts_overall_mcool": 0
            }

        # --- Setup Output Paths (UPDATED FOR HIERARCHICAL STRUCTURE) ---
        base_dir_of_mcool = os.path.dirname(mcool_full_path_abs)
        if not base_dir_of_mcool:
            base_dir_of_mcool = os.getcwd()

        results_root_folder = os.path.join(base_dir_of_mcool, "TADcaller_Results", f"TADs_{tissue}")
        contact_matrix_output_dir = os.path.join(results_root_folder, "contact_matrix")
        os.makedirs(contact_matrix_output_dir, exist_ok=True)
        print(f"Python (Debug): Contact matrix output directory: '{contact_matrix_output_dir}'")

        # --- Cooler Operations ---
        clr_path = f"{mcool_full_path_abs}::/resolutions/{resolution_bp}"
        print(f"Python (Debug): Attempting to open Cooler: {clr_path}")
        clr = cooler.Cooler(clr_path)

        # Get all chromosome names
        all_chromosomes = clr.chromnames
        print(f"Python (Debug): Found chromosomes: {all_chromosomes}")

        # Get total contacts for the *entire* cooler object once
        total_contacts_overall = int(clr.matrix(balance=False, sparse=True)[:].data.sum())


        fetch_object = clr.matrix(balance=False)

        for chrom in all_chromosomes:
            chrom_info = {"chrom": chrom}
            try:
                print(f"Python (Debug): Attempting to fetch matrix for chromosome: '{chrom}'")
                submatrix = fetch_object.fetch(chrom)

                output_path = os.path.join(contact_matrix_output_dir,
                                           f"{tissue}_{chrom}_{resolution_kb}kb_contact_matrix.txt")
                print(f"Python (Debug): Output matrix for {chrom} will be saved to: {output_path}")
                np.savetxt(output_path, submatrix, delimiter="\t", fmt="%f")

                chrom_info["output_file"] = output_path
                chrom_info["success"] = True
                chrom_info["message"] = f"✔️ Contact matrix for {chrom} saved."
            except KeyError as ke:
                overall_success = False
                chrom_info["success"] = False
                chrom_info["message"] = f"❌ Chromosome '{chrom}' not found at {resolution_kb}kb resolution. {str(ke)}"
                chrom_info["error"] = f"Chromosome not found: {str(ke)}"
            except Exception as e:
                overall_success = False
                chrom_info["success"] = False
                chrom_info["message"] = f"❌ Failed to process chromosome '{chrom}' due to an unexpected error: {str(e)}"
                chrom_info["error"] = "Python runtime error during chromosome processing"
            finally:
                processed_chromosomes_info.append(chrom_info)

        overall_message = "All requested chromosomes processed." if overall_success else "Some chromosomes failed to process."

        return {
            "success": overall_success,
            "chromosomes_processed": processed_chromosomes_info,
            "overall_message": overall_message,
            "total_contacts_overall_mcool": total_contacts_overall
        }

    except Exception as e:
        return {
            "success": False,
            "overall_message": f"❌ An unrecoverable error occurred during initial setup or Cooler object creation: {str(e)}",
            "chromosomes_processed": [],
            "total_contacts_overall_mcool": total_contacts_overall, # May be 0 if error before calculation
            "error": "Critical Python runtime error"
        }


# #-----------------------------
# ## test code for the function

# ##import os
# ##import numpy as np
# ##import cooler

# # Define the test parameters
# path_mcool_file = r"C:\Users\ghaza\Documents\ghazal\Bioinformatik_Fächer\Masterarbeit_Project\Scripts\PythonProject\test\4DNFIIH3SM5N.mcool"
# tissue = "CAKI2"
# resolution_kb = 25

# print(f"--- Starting test for {tissue} at {resolution_kb}kb resolution ---")
# print(f"Using .mcool file path: {path_mcool_file}")

# # Call the function
# result = process_contact_matrix_all_chroms(tissue, resolution_kb, path_mcool_file)

# # --- Print the results ---
# print("\n--- Function Output ---")
# print(f"Overall Success: {result.get('success')}")
# print(f"Overall Message: {result.get('overall_message')}")
# if 'total_contacts_overall_mcool' in result:
#     print(f"Total Contacts in .mcool file: {result['total_contacts_overall_mcool']}")

# print("\n--- Chromosomes Processed Details ---")
# if result.get('chromosomes_processed'):
#     for chrom_info in result['chromosomes_processed']:
#         print(f"  Chromosome: {chrom_info.get('chrom')}")
#         print(f"    Success: {chrom_info.get('success')}")
#         print(f"    Message: {chrom_info.get('message')}")
#         if chrom_info.get('output_file'):
#             print(f"    Output File: {chrom_info.get('output_file')}")
#         if chrom_info.get('error'):
#             print(f"    Error: {chrom_info.get('error')}")
#         print("-" * 30)
# else:
#     print("No chromosomes were processed or detailed information is unavailable.")

# print("\n--- End of Test ---")
# print(result)###

