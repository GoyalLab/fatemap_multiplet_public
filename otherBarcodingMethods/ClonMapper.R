#Adopting ClonMapper (GutierrezEtAl Nat Cancer 2021) data for ZhangMelzerEtAl Cell Genom 2023
#Created 20231106 by Madeline E Melzer
#Last updated 20231119 by Madeline E Melzer

library(dplyr)


# for Mac
#data <- read.delim("/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/ClonMapper/data-for-melzer/FM1.cell_record_labeled.barcode.tsv", header = FALSE)
# for PC
data <- read.delim("C:\\Users\\madel\\OneDrive - Northwestern University\\Arispe and Goyal Labs\\ZhangMelzerEtAl\\data\\ClonMapper\\data-for-melzer\\TP0_5.cell_record_labeled.barcode.tsv", header = FALSE)

# Step 2: Remove the Illumina read info column (assuming it's the first column)
data <- data[, -1] # Removes the first column

# Renaming the columns for clarity
colnames(data) <- c("UMI", "cellID", "barcode")

# Step 3: Group by cellID and UMI, select read with maximum barcode length (max = 20 nt but not all reads are complete)
result <- data %>%
  group_by(cellID, UMI) %>%
  slice_max(nchar(barcode), with_ties = FALSE) %>%
  ungroup()

# Step 4: Write the result to a new .csv file
# Replace path with the desired output file path
# Mac
#write.csv(result, "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/ZhangMelzerEtAl/data/ClonMapper/data-for-melzer/FM1_umi_counts.csv", row.names = FALSE)
# PC
write.csv(result, "C:\\Users\\madel\\OneDrive - Northwestern University\\Arispe and Goyal Labs\\ZhangMelzerEtAl\\data\\ClonMapper\\data-for-melzer\\TP0_5_umi_counts.csv", row.names = FALSE)



## combining and saving

library(readr)
library(dplyr)
library(stringr)

# Set your working directory to the folder containing the CSV files
setwd("C:\\Users\\madel\\OneDrive - Northwestern University\\Arispe and Goyal Labs\\ZhangMelzerEtAl\\data\\ClonMapper\\data-for-melzer\\")

# List all .csv files
file_list <- list.files(pattern = "\\.csv$")

# Initialize an empty list to store data frames
all_data <- list()

# Loop over the list of files
for(file_name in file_list) {
  # Read the current file
  df <- read_csv(file_name)

  # Extract the sample name from the file name using str_extract and regex
  sample_name <- str_extract(file_name, "^(.*?)_umi")

  # Add the sample column
  df$sample <- sample_name

  # Remove the UMI column
  df <- select(df, -UMI)

  # Add the data frame to the list
  all_data[[sample_name]] <- df
}

# Combine all data frames into one
final_df <- bind_rows(all_data)

# Write the final data frame to a CSV file
#write_csv(final_df, "all_samples_barcode_counts.csv")


### HERE is where I am filtering the dataframe to include only 20 nt barcodes and correct sample 

twenty_nt_barcodes <- final_df %>%
filter(nchar(barcode) == 20)

# Count the unique cellIDs associated with these 20 nt barcodes
#unique_cellIDs_with_20nt <- twenty_nt_barcodes %>%
#distinct(cellID) %>%
#nrow()  # Or use n_distinct(cellID) if you need to ensure they are unique

# Print the number of unique cellIDs
#print(unique_cellIDs_with_20nt)

modified_twenty_nt_barcodes <- twenty_nt_barcodes %>%
mutate(sample = substr(sample, 1, nchar(sample) - 4))

# View the modified dataframe
head(modified_twenty_nt_barcodes)

write_csv(modified_twenty_nt_barcodes, "all_samples_barcode_counts.csv")














#### Checking if their UMI label is actually the 10X cell ID label

#for PC
data <- read.delim("C:\\Users\\madel\\OneDrive - Northwestern University\\Arispe and Goyal Labs\\ZhangMelzerEtAl\\data\\ClonMapper\\data-for-melzer\\FM1.cell_record_labeled.barcode.tsv", header = FALSE)

# Step 2: Remove the Illumina read info column (assuming it's the first column)
data <- data[, -1] # Removes the first column

# Renaming the columns for clarity
colnames(data) <- c("cellID", "UMI", "barcode")

# reading in the list of 10X barcodes aka cellIDs

tenx = read.delim("G:\\.shortcut-targets-by-id\\1-D5WmOkOyy8I-wVx8VYZ-NDotfIYludl\\ZhangMelzerEtAl\\data\\ClonMapper\\10X\\FM1\\barcodes.tsv.gz")
tenx = data.frame(lapply(tenx, function(x) gsub("-1", "", x)), stringsAsFactors = FALSE)

# Find the intersection
intersecting_values <- intersect(data$UMI, tenx[ , 1])

# Count the number of intersecting values
number_of_intersecting_values <- length(unique(intersecting_values))

# Print the result
print(number_of_intersecting_values)





#### Old analyses

#

library(dplyr)
library(stringr)

# Assuming your data frame is named 'final_df' and the barcode column is named 'barcode'
barcode_analysis <- final_df %>%
  # Keep only 20 character barcodes
  filter(nchar(barcode) == 20) %>%
  # Add columns for the first and last 10 characters of the barcode
  mutate(
    First10 = substr(barcode, 1, 10),
    Last10 = substr(barcode, 11, 20)
  ) %>%
  # Group by the last 10 characters
  group_by(Last10) %>%
  # Summarize the number of unique first 10 characters
  summarise(
    UniqueFirst10Count = n_distinct(First10)
  ) %>%
  # Filter groups where there are multiple unique first 10 characters
  filter(UniqueFirst10Count > 1) %>%
  # Count such groups
  summarise(
    NumberOfBarcodesWithSameLast10 = n()
  )

# Print the result
print(barcode_analysis)





library(dplyr)

# Assuming your data frame is named 'final_df' and the barcode column is named 'barcode'
barcode_collapsed_info <- final_df %>%
  # Keep only 20 character barcodes
  filter(nchar(barcode) == 20) %>%
  # Add a column for the last 10 characters of the barcode
  mutate(Last10 = substr(barcode, 11, 20)) %>%
  # Group by the last 10 characters
  group_by(Last10) %>%
  # Summarize the number of distinct full barcodes for each last 10 characters
  summarise(UniqueFullBarcodes = n_distinct(barcode)) %>%
  # Filter to find groups where more than one unique full barcode maps to the same last 10 characters
  filter(UniqueFullBarcodes > 1) %>%
  # Summarise to count the total number of such full barcodes
  summarise(TotalCollapsedBarcodes = sum(UniqueFullBarcodes))
print(barcode_collapsed_info)# Print the result


# Filter to keep only 20 character barcodes
twenty_char_barcodes <- final_df %>%
  filter(nchar(barcode) == 20)

# Check the intermediate dataframe
print(head(twenty_char_barcodes))

# Add a column for the last 10 characters of the barcode
with_last_10 <- twenty_char_barcodes %>%
  mutate(Last10 = substr(barcode, 11, 20))

# Check the intermediate dataframe
print(head(with_last_10))

# Group by the last 10 characters and summarize
grouped_by_last_10 <- with_last_10 %>%
  group_by(Last10) %>%
  summarise(UniqueFullBarcodes = n_distinct(barcode))

# Check the intermediate dataframe
print(head(grouped_by_last_10))

# Filter groups with more than one unique full barcode
collisions <- grouped_by_last_10 %>%
  filter(UniqueFullBarcodes > 1)

# Check the intermediate dataframe
print(head(collisions))

# Summarise to count the total number of such full barcodes
total_collapsed_barcodes <- collisions %>%
  summarise(TotalCollapsedBarcodes = sum(UniqueFullBarcodes))

# Check the final result
print(total_collapsed_barcodes)

total_cells = barcode_data['cellID'].nunique()




### new dataframe from final_df that contains everything from final_df exccpt adding the column of the last 10 characters in each barcode
new_final_df <- final_df %>%
  mutate(Last10 = substr(barcode, nchar(barcode) - 9, nchar(barcode)))
head(new_final_df)# Display the first few rows of the new data frame


### Create a histogram of barcode length frequency
barcode_lengths <- final_df %>%
  mutate(BarcodeLength = nchar(barcode))
ggplot(barcode_lengths, aes(x = BarcodeLength)) +
  geom_histogram(binwidth = 1, color = "black", fill = "blue") +
  labs(x = "Barcode Length (nts)", y = "Frequency", title = "Histogram of Barcode Lengths") +
  theme_minimal()


##

count_instances <- result %>%
  group_by(cellID, barcode) %>%
  summarise(Count = n()) %>%
  ungroup()

# This will give you a data frame with each unique CellID and Barcode combination and their counts
print(count_instances)

count_instances_filtered <- count_instances %>%
  filter(Count > 1)

print(count_instances_filtered)




