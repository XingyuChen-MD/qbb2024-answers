# Load necessary library
library("tidyverse")

# Load data from the specified file
df <- read_tsv("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

# Extract the SUBJECT ID from SAMPID and add it as a new column
df <- df %>%
  mutate(SUBJECT = str_extract(SAMPID, "[^-]+-[^-]+"), .before = 1)

# Count the number of samples per SUBJECT and display in descending order
df %>%
  group_by(SUBJECT) %>%
  summarise(my_counts = n()) %>%
  arrange(desc(my_counts))
# Note: K-562 has the most samples, GTEX-OHPN and GTEX-QDT8 have the least samples

# Count the number of samples per tissue type (SMTSD) and display in descending order
df %>%
  group_by(SMTSD) %>%
  summarise(sample_count = n()) %>%
  arrange(desc(sample_count))
# Note: Whole Blood has the most samples, Skin - Not Sun Exposed (Suprapubic) has the least samples

# Filter for samples from a specific subject (GTEX-NPJ8) and save to a new data frame
df_npj8 <- df %>%
  filter(SUBJECT == "GTEX-NPJ8")

# Count the number of samples per tissue type for the selected subject
df_npj8 %>%
  group_by(SMTSD) %>%
  summarise(sample_count = n()) %>%
  arrange(desc(sample_count))
# Note: Whole Blood has the most samples

# View columns 15 to 20 for the Whole Blood samples from GTEX-NPJ8
blood <- df_npj8 %>%
  filter(SMTSD == "Whole Blood")
blood[, 15:20]

# Filter out rows where SMATSSCR is not NA
df_na <- df %>%
  filter(!is.na(SMATSSCR))

# Calculate the mean SMATSSCR score for each subject and filter those with a mean score of 0
df_0 <- df_na %>%
  group_by(SUBJECT) %>%
  summarise(mean_SMATSSCR = mean(SMATSSCR)) %>%
  filter(mean_SMATSSCR == 0)
# Note: 15 subjects have a mean SMATSSCR score of 0

# Calculate the range and mean of SMATSSCR mean scores across all subjects
df_mean <- df_na %>%
  group_by(SUBJECT) %>%
  summarise(mean_SMATSSCR = mean(SMATSSCR))

range(df_mean$mean_SMATSSCR) # Range of SMATSSCR means: [0, 2.571429]
mean(df_mean$mean_SMATSSCR)  # Mean of SMATSSCR means: 0.9050125

# Visualize the distribution of SMATSSCR mean scores using a dot plot
ggplot(df_mean, aes(x = mean_SMATSSCR)) +
  geom_dotplot(binwidth = 0.1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(x = "Mean SMATSSCR Score", y = "Frequency", title = "Distribution of Mean SMATSSCR Scores Across Subjects") +
  theme_minimal(base_size = 12) +
  theme(
    title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
