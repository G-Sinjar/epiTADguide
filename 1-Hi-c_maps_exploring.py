import os
import numpy as np
import cooler


def process_contact_matrix(tissue, chrom, resolution_kb, path_mcool_file):
    """
    Args:
        tissue (str): The tissue name (e.g., "CAKI2").
        chrom (str or int): The chromosome name (e.g., "chr2", "chrX").
        resolution_kb (int): The resolution in kilobases (e.g., 25).
        path_mcool_file (str): The path to the .mcool file. Can be:
                                - Full path with .mcool extension (e.g., "/data/file.mcool")
                                - Path without .mcool extension (e.g., "/data/file")
                                - Directory (will raise an error)
    Returns:
        dict: A dictionary containing:
              - "success" (bool): True if successful, False otherwise.
              - "output_file" (str, optional): Path to the saved contact matrix file.
              - "message" (str): A status or error message.
              - "error" (str, optional): A brief error description if not successful.
              - "total_contacts" (int, optional): Sum of all contacts in the *entire* .mcool file (not just the fetched chromosome).
    """
    try:
        # --- Input Normalization and Validation ---
        # Ensure resolution_kb is an integer
        try:
            resolution_kb = int(resolution_kb)
            resolution_bp = resolution_kb * 1000
        except ValueError:
            return {
                "success": False,
                "message": f"❌ Invalid resolution input: '{resolution_kb}'. Resolution must be an integer in kilobases (e.g., 25).",
                "error": "Invalid resolution input"
            }

        # Determine the full .mcool file path and check for existence
        mcool_full_path = path_mcool_file
        if not mcool_full_path.endswith(".mcool"):
            mcool_full_path = f"{path_mcool_file}.mcool"

        # Resolve to an absolute path for robust checking
        mcool_full_path_abs = os.path.abspath(mcool_full_path)

        if not os.path.exists(mcool_full_path_abs):
            return {
                "success": False,
                "message": f"❌ File not found: '{mcool_full_path_abs}'. Please check the path and filename (with or without .mcool extension) and try again.",
                "error": "Missing .mcool file"
            }

        # Critical check: Ensure it's a file, not a directory
        if os.path.isdir(mcool_full_path_abs):
            return {
                "success": False,
                "message": f"❌ Input path resolves to a directory: '{mcool_full_path_abs}'. \nPlease provide the full path to the .mcool file, including its name (e.g., /data/my_file.mcool).",
                "error": "Directory instead of file path"
            }

        # --- Setup Output Paths (UPDATED FOR HIERARCHICAL STRUCTURE) ---
        # The base directory for results will be the directory where the .mcool file is located
        base_dir_of_mcool = os.path.dirname(mcool_full_path_abs)
        if not base_dir_of_mcool:  # If empty (file is in CWD)
            base_dir_of_mcool = os.getcwd()

        # Define the root for all results for this run: TADcaller_Results/TADs_YourTissue/
        # This aligns with the R Shiny app's definition for consistent folder creation
        results_root_folder = os.path.join(base_dir_of_mcool, "TADcaller_Results", f"TADs_{tissue}")

        # Create the contact_matrix subfolder within this root
        contact_matrix_output_dir = os.path.join(results_root_folder, "contact_matrix")
        os.makedirs(contact_matrix_output_dir, exist_ok=True)
        print(f"Python (Debug): Contact matrix output directory: '{contact_matrix_output_dir}'")


        # --- Cooler Operations ---
        clr_path = f"{mcool_full_path_abs}::/resolutions/{resolution_bp}"
        print(f"Python (Debug): Attempting to open Cooler: {clr_path}")
        # Instantiate cooler object
        clr = cooler.Cooler(clr_path)

        # Extract and save contact matrix
        fetch_object = clr.matrix(balance=False)

        print(f"Python (Debug): Attempting to fetch matrix for chromosome: '{chrom}'")
        # Attempt to fetch the submatrix - this is where KeyError might occur
        submatrix = fetch_object.fetch(chrom)

        # Output path for the contact matrix (uses the new contact_matrix_output_dir)
        output_path = os.path.join(contact_matrix_output_dir,
                                   f"{tissue}_{chrom}_{resolution_kb}kb_contact_matrix.txt")
        print(f"Python (Debug): Output matrix will be saved to: {output_path}")
        np.savetxt(output_path, submatrix, delimiter="\t", fmt="%f")

        # Get contact count for the *entire* cooler object
        total_contacts = clr.matrix(balance=False, sparse=True)[:].data.sum()

        return {
            "success": True,
            "output_file": output_path, # This is the path to the contact matrix
            "message": f"✔️ Contact matrix for {chrom} at {resolution_kb}kb resolution saved to: '{output_path}'.",
            "total_contacts": int(total_contacts)
        }

    except KeyError as ke:
        # This specific error usually means the chromosome was not found in the .mcool file
        return {
            "success": False,
            "message": f"❌ Chromosome '{chrom}' was not found in the .mcool file at {resolution_kb}kb resolution. Please check the chromosome name.",
            "error": f"Chromosome not found: {str(ke)}"
        }
    except Exception as e:
        # Catch any other unexpected errors
        return {
            "success": False,
            "message": f"❌ Failed to process contact matrix due to an unexpected error: {str(e)}",
            "error": "Python runtime error"
        }


    
# ## trying the function
# result = process_contact_matrix(
#     tissue="CAKI2",
#     chrom="chrX",
#     resolution_kb = 25,
#     path_mcool_file="./test/4DNFIIH3SM5N"
# )
#
# if result["success"]:
#     print(result["message"])
#     print(f"Total contacts: {result['total_contacts']:,}")
# else:
#     print(result["message"])
#     print("Error:", result["error"])