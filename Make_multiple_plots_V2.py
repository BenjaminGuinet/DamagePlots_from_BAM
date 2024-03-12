import argparse
import os
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import pysam
import re,os
from tabulate import tabulate
import PyPDF2
import numpy as np 
import subprocess
import matplotlib.ticker as ticker
import glob
import statistics
import pysam
import os

# module load picard/3.1.1
# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Split bam files according to contigs for each species     .\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-b", "--bam",help="The bam file name")
args = parser.parse_args()


# Usage eg : python3 Make_multiple_plots_V2.py -b sorted_Mammuthus-MD002_TH1177_sp014523485_RdrVT6VMmm_extracted.bam 
UDG_treated = "yes"

bam_file_path = args.bam

# Check if BAM file has an index

def check_and_index_bam(bam_path):
    # The index file is usually the same name as the BAM file but with .bai extension
    index_file_path = bam_path + '.bai'
    # Check if the index file exists
    if os.path.exists(index_file_path):
        print("Index file already exists.")
    else:
        print("Index file does not exist. Creating index...")
        # Create the BAM file index. This function automatically saves the index file
        # in the same directory as the BAM file.
        pysam.index(bam_path)
        print("Index created and saved to", index_file_path)

check_and_index_bam(bam_file_path)

def get_read_length_distribution(bam_file):
    read_lengths = {}
    read_lengths_all = []
    if not os.path.exists(bam_file+".bai"):
        command = ["samtools index ", bam_file]
        subprocess.run(" ".join(command), shell=True, check=True)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            length = read.query_length
            read_lengths_all.append(length)
            if length in read_lengths:
                read_lengths[length] += 1
            else:
                read_lengths[length] = 1
    # If the proportion of reads length <100 is > 90% of the total number of reads, score = 1
    if sum([read_lengths[length] for length in read_lengths if length < 100]) / sum(read_lengths.values()) > 0.9:
        score = 1
    else:
        score = 0
    mean_read_length = sum(read_lengths_all) / len(read_lengths_all)
    # Get standard deviation of read lengths
    if len(read_lengths_all) > 1:  # Standard deviation calculation requires at least 2 values
        stdev_read_length = round(statistics.stdev(read_lengths_all), 2)
    else:
        stdev_read_length = 0.0
    # Return a tuple of the dictionary and the mean read length
    return read_lengths, mean_read_length, score, stdev_read_length

def get_identity_distribution_and_avg_mapq(bam_file):
    identity_values = []
    total_nb_reads = 0
    total_mapq = []  # To store the sum of mapping qualities
    if not os.path.exists(bam_file + ".bai"):
        command = ["samtools index", bam_file]
        subprocess.run(" ".join(command), shell=True, check=True)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not read.is_unmapped:
                identity = read.get_tag("NM")  # NM tag represents edit distance
                read_length = read.query_length
                percentage_identity = ((read_length - identity) / read_length) * 100
                identity_values.append(percentage_identity)
                total_mapq.append(read.mapping_quality)  
                total_nb_reads += 1
    avg_mapq = sum(total_mapq) / len(total_mapq)
    avg_mapq_rounded = round(avg_mapq, 2)
    # Get mean of identity values
    avg_identity = sum(identity_values) / len(identity_values)
    # Calculate standard deviation for identity values
    if len(identity_values) > 1:  # Standard deviation calculation requires at least 2 values
        stdev_identity = round(statistics.stdev(identity_values), 2)
    else:
        stdev_identity = 0.0  # Standard deviation is not defined for a single value
    # Calculate the standard deviation of the mapping quality
    if len(total_mapq) > 1:  # Standard deviation calculation requires at least 2 values
        stdev_mapq = round(statistics.stdev(total_mapq), 2)
    else:
        stdev_mapq = 0.0  # Standard deviation is not defined for a single value
    return identity_values, total_nb_reads, avg_mapq_rounded,stdev_identity, stdev_mapq, avg_identity

def plot_read_length_distribution(read_length_distribution, mean_read_length, PMD_read_length_distribution, mean_PMD_read_length_distribution):
    plt.bar(read_length_distribution.keys(), read_length_distribution.values(), width=1.0, color='skyblue', label='Reads')
    plt.bar(PMD_read_length_distribution.keys(), PMD_read_length_distribution.values(), width=1.0, color='red', alpha=0.5, label='Reads PMDscore > 3')
    # Plotting red dashed line for the mean of normal reads
    plt.axvline(x=mean_read_length, color='black', linestyle='--', linewidth=1.5)
    # Displaying the mean value of normal reads as text
    plt.text(mean_read_length + 0.5, plt.ylim()[1] * 0.9, f"Mean : {mean_read_length:.2f}", color='black')
    # Plotting blue dashed line for the mean of PMD reads
    plt.axvline(x=mean_PMD_read_length_distribution, color='red', linestyle='--', linewidth=1.5)
    # Displaying the mean value of PMD reads as text
    plt.text(mean_PMD_read_length_distribution + 0.5, plt.ylim()[1] * 0.85, f"Mean : {mean_PMD_read_length_distribution:.2f}", color='red')
    plt.xlabel('Read Length')
    plt.ylabel('Count')
    plt.title('Read Length Distribution')
    plt.legend()

def plot_identity_distribution(identity_values):
    plt.hist(identity_values, bins=30, color='lightgreen', edgecolor='black')
    # Calculate mean of identity values
    mean_identity = np.mean(identity_values)
    # Plotting red dashed line for the mean
    plt.axvline(x=mean_identity, color='black', linestyle='--', linewidth=1.5)
    # Displaying the mean value as text
    plt.text(mean_identity + 0.5, plt.ylim()[1] * 0.9, f"Mean: {mean_identity:.2f}", color='black')
    plt.xlabel('Percentage of Identity')
    plt.ylabel('Count')
    plt.title('Percentage of Identity Distribution')

def amber_plot(bam_file_path):
    # Run Amber 
    command = [
        "python ",
         "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/TOOLS/Amber_last.py --bamfiles=",
        bam_file_path, " --output=", re.sub(".bam", "_amber", bam_file_path)]
        # Run the command using subprocess
    subprocess.run(''.join(command), shell=True)
     ### Open amber output 
    # Path to your text file
    file_path = 'sortedP033_nonUDG_treated_Mesorhizobium_manganicum_WmbN6dKam7_extracted_amber.txt'
    # Read the entire file into a single string
    with open(re.sub(".bam", "_amber", bam_file_path)+".txt", 'r') as file:
            file_contents = file.read()
        # Identify the start and end of the section of interest
    start_marker = "DISTANCE_FROM_READ_END\tC-to-T\tCpG-to-TG\tother"
    end_marker = "----------------------------------------------------------------------------------------------------"
    start_index = file_contents.find(start_marker)
    end_index = file_contents.find(end_marker, start_index)
    # Extract the section of interest
    section_of_interest = file_contents[start_index:end_index].strip()
    # Convert the section to a DataFrame
    data = [line.split('\t') for line in section_of_interest.split('\n')]
    tab = pd.DataFrame(data[1:], columns=data[0])
    # Convert numeric columns from string to float
    numeric_cols = tab.columns[1:]  # Assuming the first column is non-numeric
    tab[numeric_cols] = tab[numeric_cols].apply(pd.to_numeric, errors='coerce')
    # Rename columns
    tab.columns = ["distance_from_read_end", "C-to-T", "CpG-to-TpG", "other"]
    if tab["CpG-to-TpG"][0] > 0.05:
        score = 1
    else:
        score = 0
    # USe pydamage to compute a pvalue 
    print("Running pydamage...")
    pydamage_dir="pydamage_"+re.sub(".bam","",os.path.basename(bam_file_path))
    cmd = [
    "/home/benjguin/.conda/envs/My_conda/bin/pydamage  --outdir ",pydamage_dir,  # Path to the pydamage executable
     " analyze ",  # Command to analyze
    bam_file_path,  # Input file
    " --group "]  # Group argument
    subprocess.run(''.join(cmd), shell=True)
    print("pydamage done.")
    # Extracting data from the DataFrame
    # Extracting data from the DataFrame
    if UDG_treated == "yes":
        y1 = tab["C-to-T"]
        y2 = tab["CpG-to-TpG"]
        y3 = tab["other"]
    else:
        y1 = tab["CpG-to-TpG"]
        y2 = tab["C-to-T"]
        y3 = tab["other"]
    x1 = tab["distance_from_read_end"]
    x2 = tab["distance_from_read_end"]
    x3 = tab["distance_from_read_end"]
    # Extract title from file name
    # Create the DNA damage subplot within the larger subplot grid
    plt.xticks(list(range(0, 31, 2)))  # Replace ax.set_xticks with plt.xticks
    all_y = pd.concat([y1, y2, y3])
    ymin, ymax = all_y.min(), all_y.max()
    ymin -= (ymax - ymin) * 0.1  # Optional: decrease ymin by 10% of the range for padding
    ymax += (ymax - ymin) * 0.1  # Optional: increase ymax by 10% of the range for padding
    ymax += 0.10
    plt.ylim(ymin, ymax)  # Replace ax.set_ylim with plt.ylim
    plt.xlabel("distance from read end (bp)")  # Replace ax.set_xlabel with plt.xlabel
    plt.ylabel("mismatch frequency")  # Replace ax.set_ylabel with plt.ylabel
    # Put the line y1 in grey and y2 in red
    if UDG_treated == "yes":
        plt.plot(x1, y1, label="C to T", linewidth=1.5, linestyle='dashed',color="black")  # Replace ax.plot with plt.plot
        plt.plot(x2, y2, label="CpG to TpG", linewidth=1.5, color="red")  # Replace ax.plot with plt.plot
        plt.plot(x3, y3, label="other", linestyle='dotted',linewidth=1.5, color="grey")  # Replace ax.plot with plt.plot
    else:
        plt.plot(x1, y1, label="CpG to TpG", linewidth=1.5, linestyle='dashed',color="grey")  # Replace ax.plot with plt.plot
        plt.plot(x2, y2, label="C to T", linewidth=1.5, color="red")  # Replace ax.plot with plt.plot
        plt.plot(x3, y3, label="other", linestyle='dotted',linewidth=1.5, color="grey")  # Replace ax.plot with plt.plot
    # Add title
    plt.title("DNA damage pattern")  # Replace ax.set_title with plt.title
    plt.legend(loc="upper right", edgecolor="black")  # No change needed here
    return plt, score, tab

def evenness_plot(bam_file_path):
    depth_data = list(pysam.depth('-a', bam_file_path, split_lines=True))
    split_data = [entry.split('\t') for entry in depth_data]
    # Convert the list of lists into a DataFrame
    df = pd.DataFrame(split_data, columns=['Reference', 'Position', 'Depth'])
    # add Position and Depth as numeric values
    df['Position'] = pd.to_numeric(df['Position'])
    df['Depth'] = pd.to_numeric(df['Depth'])
    Mean_coverage = round(df["Depth"].mean(), 2)
    # Calculate the percentage of base with cov > 0 
    Percentage_bases_covered = round((df["Depth"] > 0).sum() / len(df["Depth"]) * 100, 2)
    N_tiles = 100
    min_val = df["Position"].min(numeric_only=True)
    max_val = df["Position"].max(numeric_only=True)
    step = (max_val - min_val) / N_tiles
    tiles = np.arange(0, N_tiles + 1) * step
    V4 = np.array([])
    V5 = np.array([])
    Nb_window  =0
    for i in range(len(tiles) - 1):
        Nb_window +=1
        df_temp = df[(df["Position"] >= tiles[i]) & (df["Position"] < tiles[i + 1])]
        temp_V4 = np.repeat((df_temp["Depth"] > 0).sum() / len(df_temp["Depth"]), len(df_temp["Depth"]))
        temp_Nb_window = np.repeat(Nb_window, len(df_temp["Depth"]))
        temp_V4[np.isnan(temp_V4)] = 0  # Replace NaN in temp_V4 with 0
        V4 = np.append(V4, temp_V4)
        V5 = np.append(V5, temp_Nb_window)
    # Create a new array by repeating the last element twice
    try:
        df["V4"] = V4
        df["V5"] = V5
    except:
     try:
        V4= np.append(V4, V4[-1])
        V5= np.append(V5, V5[-1])
        df["V4"] = V4
        df["V5"] = V5
     except:
      try:
        diff=abs(len(df)-len(V4))
        for i in range(1,diff+1):
          V4= np.append(V4, V4[-1])
          V5= np.append(V5, V5[-1])
        df["V4"] = V4
        df["V5"] = V5
      except:
        diff=abs(len(df)-len(V4))
        V4= V4[:-diff]
        V5 = V5[:-diff]
        df["V4"] = V4
        df["V5"] = V5
    # remove duplicated values within columns V5 
    df_sub= df.drop_duplicates(subset=['V5'], keep='first')
    # Count number of V4 values < 0.01
    Percentage_window_covered= round((1-np.sum(df_sub["V4"] < 0.01) / len(df_sub))*100, 2)
    # Split reference genome into tiles, compute breadth of coverage for each tile
    N_tiles = 500
    step = (df['Position'].max() - df['Position'].min()) / N_tiles
    tiles = np.arange(0, N_tiles + 1) * step
    boc = []
    Nb_window  =0
    V6 = np.array([])
    for i in range(len(tiles) - 1):
        Nb_window +=1
        df_loc = df[(df['Position'] >= tiles[i]) & (df['Position'] < tiles[i+1])]
        temp_Nb_window = np.repeat(Nb_window, len(df_loc["Depth"]))
        boc.extend(np.repeat(np.sum(df_loc['Depth'] > 0) / len(df_loc['Depth']), len(df_loc)))
        V6 = np.append(V6, temp_Nb_window)
    try:
        df["V6"] = V6
    except:
     try:
        V6= np.append(V6, V6[-1])
        df["V6"] = V6
     except:
      try:
        diff=abs(len(df)-len(V6))
        for i in range(1,diff+1):
          V6= np.append(V6, V6[-1])
        df["V6"] = V6
      except:
        diff=abs(len(df)-len(V6))
        V6= V6[:-diff]
        df["V6"] = V6
    # Handling NaN values in 'boc'
    boc = [0 if np.isnan(x) else x for x in boc]
    # Add two times the last value of boc in the boc list
    try:
     df['boc'] = boc
    except:
     try:
      boc.extend([boc[-1]] * 1)
      df['boc'] = boc
     except:
      try:
       boc.extend([boc[-1]] * 1)
       df['boc'] = boc
      except:
       boc.extend([boc[-1]] * 1)
       df['boc'] = boc
    # Plotting
    # put  df['Position2']  as the row number 
    df['Position2'] = df.index
    grouped = df.groupby('V6')
    for name, group in grouped:
        plt.plot(group['Position2'], group['boc'], linestyle='-', linewidth=1, color='black', label='Lines')
    # Linking lines between last of one group to the first of the next
    for i in range(1, len(grouped)):
        prev_group = grouped.get_group(i)
        next_group = grouped.get_group(i + 1)
        plt.plot([prev_group.iloc[-1]['Position2'], next_group.iloc[0]['Position2']],
                [prev_group.iloc[-1]['boc'], next_group.iloc[0]['boc']],
                color='gray', linestyle='-',linewidth=0.5)
    plt.xlabel('Reference genome position')
    plt.ylabel('Fraction of covered genome')
    # Add in the title the percentage of windows covered
    plt.title(f"{Percentage_window_covered}% windows & {Percentage_bases_covered}% genome covered")  # Adjust title font siz
    plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: '{:0.0e}'.format(x)))
    # if Percentage_window_covered > 50, add score of 1, if > 95, add score of 2
    if Percentage_window_covered > 95:
        score = 2
    elif Percentage_window_covered > 50:
        score = 1
    else:
        score = 0
    return plt, Mean_coverage, score, Percentage_window_covered, Percentage_bases_covered
    ### 

def PMDscore_run(bam_file_path):
    extracted_bam_file_PMDscores_extracted = re.sub("\.bam", "_PMDscores_extracted_sup3_CpG.bam", bam_file_path)
    command1 = ["samtools view -h ", bam_file_path," | python2 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/TOOLS/pmdtools.0.60.py --CpG  --threshold 3 --header | samtools view -Sb - ", ">", extracted_bam_file_PMDscores_extracted]
    subprocess.run(" ".join(command1), shell=True, check=True)
    command1 = ["samtools index ", extracted_bam_file_PMDscores_extracted]
    subprocess.run(" ".join(command1), shell=True, check=True)
    get_read_length_distribution_output_PMD=get_read_length_distribution(extracted_bam_file_PMDscores_extracted)
    # Run to get the distribution of PMD scores 
    extracted_bam_file_PMDscores_tab = re.sub("\.bam", "_PMDscores_CpG.txt", bam_file_path)
    command1 = ["samtools", "view -h", bam_file_path, "| python2 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/TOOLS/pmdtools.0.60.py --CpG --printDS ", ">", extracted_bam_file_PMDscores_tab]
    subprocess.run(" ".join(command1), shell=True, check=True)
    return get_read_length_distribution_output_PMD, extracted_bam_file_PMDscores_extracted,extracted_bam_file_PMDscores_tab

def PMDscore_plot(tab_bam_file_PMDscores):
    pmd_scores = pd.read_csv(tab_bam_file_PMDscores, header=None, sep="\t")
    Percentage_score_sup_3 = round(np.sum(pmd_scores[3] > 3) / len(pmd_scores) * 100, 2)
    Percentage_score_sup_1pt5 = round(np.sum(pmd_scores[3] > 1.5) / len(pmd_scores) * 100, 2)
    # If the proportion of PMD scores > 3 is > 10% of the total number of reads, score = 1
    if Percentage_score_sup_3 > 10:
        score = 1
    else:
        score = 0
    # Create a histogram
    # remove values == 0 in the list pmd_scores[3]
    pmd_scores_values  = pmd_scores[3][pmd_scores[3] != 0]
    plt.hist(pmd_scores_values , bins=500, color="black")
    plt.xlabel("PMDscores")
    plt.ylabel("Frequency")
    plt.title("Histogram of PMD Scores")
    # Set Frequenxy max as 150 
    plt.xlim(min(pmd_scores[3]), 15)
    # Adding text for Percentage_score_sup_3 in the top right corner
    text = f"PMDscore CpG > 3: {Percentage_score_sup_3}%"
    plt.text(0.6, 0.9, text, ha='center', va='center', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.8))
    # Add a vertical line at x=3
    plt.axvline(x=3, color='red', linestyle='--', linewidth=1.5)
    return plt, score, Percentage_score_sup_3, Percentage_score_sup_1pt5

for i in ["A"]:
            print("Creating plots for the bamfile : ",bam_file_path)
            Name=bam_file_path
            # First wee need to process the bam file to remove duplicates and 
            # Only select reads with maping quality > 0 
            cmd = ["samtools view -h -q 1 ", bam_file_path, " > ", "quality_filtred_",bam_file_path]  # Group argument
            subprocess.run(''.join(cmd), shell=True)
            # Remove all bad headers 
            cmd1 = ["samtools sort ","quality_filtred_",bam_file_path, "> sorted_","quality_filtred_",bam_file_path ]  
            subprocess.run(''.join(cmd1), shell=True)
            cmd2 = ["/crex/proj/snic2022-6-144/nobackup/BENJAMIN/my_conda/bin/filterBAM --bam ","sorted_quality_filtred_",bam_file_path, " --min-read-length 30 --min-read-ani 0 --min-avg-read-ani 0 --min-normalized-gini 1  -t 1 --stats ", "quality_filtred_",bam_file_path,"_filterBAM_stats.csv --bam-filtered ", "filterBAM_quality_filtred_",bam_file_path, " --stats-filtered ", "filterBAM_quality_filtred_",bam_file_path,"_stats.csv"]  
            subprocess.run(''.join(cmd2), shell=True)
            #cmd3 = ["samtools sort ","quality_filtred_",bam_file_path, "> sorted_","quality_filtred_",bam_file_path ]  
            #subprocess.run(''.join(cmd3), shell=True)
            cmd4 = ["samtools sort", " filterBAM_quality_filtred_",bam_file_path , " > ", "sorted_filterBAM_quality_filtred_",bam_file_path]  # Group argument
            subprocess.run(''.join(cmd4), shell=True)
            # Use subprocess to run the command to change the name of the file
            Full_name = "sorted_filterBAM_quality_filtred_"+bam_file_path
            # remove "_filtred_sorted_" from the name
            #Full_name=re.sub("_filtred_sorted_", "_", Full_name)
            cmd4 = ["java -jar $PICARD MarkDuplicates I=",Full_name, " O=NoDup_",Full_name," REMOVE_DUPLICATES=true, CREATE_INDEX=true, M=NoDup_",Full_name,".metrics"]  # Group argument
            subprocess.run(''.join(cmd4), shell=True)
            #cmd5 = ["mv ", "sorted_filterBAM_quality_filtred_",bam_file_path, " ", Full_name]  # Group argument
            #subprocess.run(''.join(cmd5), shell=True)
            #cmd6 = ["samtools index ",Full_name]  # Group argument
            #subprocess.run(''.join(cmd6), shell=True)
            # Remove all the intermediate files instead of the final bam file
            cmd7 = ["rm ", "quality_filtred_",bam_file_path]  # Group argument
            subprocess.run(''.join(cmd7), shell=True)
            cmd8 = ["rm ", "sorted_quality_filtred_",bam_file_path]  # Group argument
            subprocess.run(''.join(cmd8), shell=True)
            cmd9 = ["rm ", "filterBAM_quality_filtred_",bam_file_path]  # Group argument
            subprocess.run(''.join(cmd9), shell=True)
            cmd10 = ["rm ", "quality_filtred_",bam_file_path]  # Group argument
            subprocess.run(''.join(cmd10), shell=True)
            cmd11 = ["rm ", "filterBAM_quality_filtred_",bam_file_path,".bai"]  # Group argument
            subprocess.run(''.join(cmd11), shell=True)
            cmd12 = ["rm ", "sorted_filterBAM_quality_filtred_",bam_file_path]  # Group argument
            subprocess.run(''.join(cmd12), shell=True)
            Full_name ="NoDup_"+Full_name
            bam_file_path = Full_name
            # Sort and index the bam file 
            # Get read length distribution
            get_read_length_distribution_output = get_read_length_distribution(bam_file_path)
            identity_distribution = get_identity_distribution_and_avg_mapq(bam_file_path)
            read_length_distribution = get_read_length_distribution_output[0]
            mean_read_length= get_read_length_distribution_output[1]
            sd_read_length = get_read_length_distribution_output[3]
            print("\n")
            print("Running PMDtools...")
            get_read_length_distribution_output_PMD=PMDscore_run(bam_file_path)
            PMD_read_length_distribution=get_read_length_distribution_output_PMD[0][0]
            mean_PMD_read_length_distribution=get_read_length_distribution_output_PMD[0][1]
            extracted_bam_file_PMDscores=get_read_length_distribution_output_PMD[1]
            tab_bam_file_PMDscores=get_read_length_distribution_output_PMD[2]
            read_length_distribution_score=get_read_length_distribution_output[2]
            print("PMDtools done.")
            print("\n")
            # Get the text for the general title
            general_title_text = re.sub("_extracted", "", Name)
            plt.figure(figsize=(12, 12))
            plt.suptitle(general_title_text, fontsize=14, fontweight='bold', ha='center')
            plt.subplot(3, 2, 1)
            plot_read_length_distribution(read_length_distribution,mean_read_length, PMD_read_length_distribution, mean_PMD_read_length_distribution)
            print("plot_read_length_distribution done.")
            print("\n")
            # Replace 'identity_distribution' with your identity distribution data or function
            plt.subplot(3, 2, 2)
            plot_identity_distribution(identity_distribution[0])
            print("plot_identity_distribution done.")
            print("\n")
            plt.subplot(3, 2, 3)
            amber_plot_output=amber_plot(bam_file_path)
            # This is the plot
            amber_plot_output[0]
            # This is the score
            amber_score=amber_plot_output[1]
            amber_table=amber_plot_output[2]
            #
            del amber_table['other']
            amber_table = amber_table.melt(id_vars='distance_from_read_end', var_name='variable', value_name='value')
            # Creating new column names
            amber_table['new_column'] = amber_table['variable'] + '_' + amber_table['distance_from_read_end'].astype(str)
            # Creating the final DataFrame
            amber_table = amber_table.pivot_table(index=None, columns='new_column', values='value', aggfunc='first')
            amber_table.reset_index(drop=True, inplace=True)
            #
            print("plot_amber done.")
            print("\n")
            plt.subplot(3, 2, 4)
            evenness_plot_output= evenness_plot(bam_file_path)
            # This is the plot
            evenness_plot_output[0]
            # This is the mean coverage
            Mean_coverage=evenness_plot_output[1]
            evenness_score=evenness_plot_output[2]
            Percentage_window_covered = evenness_plot_output[3]
            Percentage_bases_covered = evenness_plot_output[4]
            print("plot_evenness_coverage done.")
            print("\n")
            plt.subplot(3, 2, 5)
            PMDcore_plot=PMDscore_plot(tab_bam_file_PMDscores)
            # This is the plot
            PMDcore_plot[0]
            PMD_score=PMDcore_plot[1]
            PMD_perc_score_sup_3=PMDcore_plot[2]
            PMD_perc_score_sup_1pt5=PMDcore_plot[3]
            print("plot_PMDscore done.")
            print("\n")
            # Get the summ of all scores
            Sum_scores=amber_score+read_length_distribution_score+evenness_score+PMD_score
            Sum_scores=str(Sum_scores)+"/5"
            # Create data for subplots
            # Create a DataFrame from the Kraken report
            N_tot_reads=identity_distribution[1]
            Average_mapq=identity_distribution[2]
            sdAverage_mapq=identity_distribution[4]
            sdAverage_ID=identity_distribution[3]
            Mean_identity=identity_distribution[5]
            pydamage_dir="pydamage_"+re.sub(".bam","",os.path.basename(bam_file_path))
            pydamage_results=pd.read_csv(pydamage_dir+"/"+"pydamage_results.csv")
            data = {"Nb reads": [N_tot_reads],
            'Mean coverage': [Mean_coverage],
            "Average mapq": str(Average_mapq) + " (sd=" + str(sdAverage_mapq) + ")",
            "Pydamage Pmax": str(pydamage_results["damage_model_pmax"].iloc[0]) + " (sd=" + str(pydamage_results["damage_model_pmax_stdev"].iloc[0]) + ")",
            "Pydamage Pvalue": str(pydamage_results["pvalue"].iloc[0]),
            "Pydamage RMSE": str(pydamage_results["RMSE"].iloc[0])}
            df = pd.DataFrame(data)
            df = df.T.reset_index()
            df.columns = ['Attributes', 'Values']
            plt.subplot(3, 2, 6)
            table = plt.table(cellText=df.values, colLabels=df.columns, loc='center', cellLoc='center')
            # Adjust table properties
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            table.scale(1.1, 1.7)  # Scale the table size
            cell_width = 0.15  # Adjust cell width
            table.auto_set_column_width([i for i in range(4)])
            table.set_fontsize(9)
            # Hide axes for the subplot with the table
            plt.gca().axis('off')
            plt.tight_layout(rect=[0, 0, 1, 0.95]) 
            print("Saving the plot to : "+bam_file_path + "_dna_damage_plots.pdf")
            plt.savefig(bam_file_path + "_dna_damage_plots.pdf")
            print("\n")
            pydamage_dir="pydamage_"+re.sub(".bam","",os.path.basename(bam_file_path))  
            pydamage_results=pd.read_csv(pydamage_dir+"/"+"pydamage_results.csv")
            # CSV OUTPUT
            csv = { 
            "Bam_file_name" : [bam_file_path],
            "General_ancient_score": [Sum_scores],
            "Nb_reads_mapped": [N_tot_reads],
            'Mean_coverage': [Mean_coverage],
            'Percentage_window_covered': [Percentage_window_covered],
            "Average_mapq": Average_mapq,
            "sdAverage_mapq": sdAverage_mapq,
            "Pydamage_Pmax": pydamage_results["damage_model_pmax"].iloc[0],
            "Pydamage_Pvalue": pydamage_results["pvalue"].iloc[0],
            "Pydamage_RMSE": pydamage_results["RMSE"].iloc[0],
            "Percentage_PMD_score_sup_1pt5" : PMD_perc_score_sup_1pt5,
            "Percentage_PMD_score_sup_3" : PMD_perc_score_sup_3,
            "sd_identity": sdAverage_ID,
            "Mean_identity": Mean_identity,
            "Mean_read_length": mean_read_length,
            "sd_read_length": sd_read_length,
            "Percentage_bases_covered": Percentage_bases_covered
            }
            # Creating a DataFrame
            output_csv = pd.DataFrame(csv)
            # Add the columns from amber_table to the DataFrame
            output_csv = pd.concat([output_csv, amber_table], axis=1)
            # Save the table to a CSV file in the specified directory
            file_path = re.sub(".bam", "_dna_damage.csv", bam_file_path)  # Adjust the file name as needed
            output_csv.to_csv(file_path, index=False)
            print ("The csv file has been saved to : "+file_path)
            print("\n")
            print("All done.")
