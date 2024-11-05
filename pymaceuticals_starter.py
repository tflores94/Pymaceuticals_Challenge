#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# - Add your analysis here.
#  

# In[161]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "Mouse_metadata.csv"
study_results_path = "Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single DataFrame
study_data_merged = pd.merge(study_results, mouse_metadata, how="left", on="Mouse ID")

# Display the data table for preview
study_data_merged.head()


# In[163]:


# Checking the number of mice.
len(study_data_merged["Mouse ID"].unique())


# In[165]:


# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint.
find_duplicated_rows = study_data_merged.duplicated(subset=['Mouse ID', 'Timepoint'])
duplicated_mice = study_data_merged.loc[find_duplicated_rows,'Mouse ID'].unique()
duplicated_mice


# In[167]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_data = study_data_merged[study_data_merged['Mouse ID'].isin(duplicated_mice)==False]
clean_data.head()


# In[169]:


# Checking the number of mice in the clean DataFrame.
len(clean_data["Mouse ID"].unique())


# ## Summary Statistics

# In[172]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen:
# mean, median, variance, standard deviation, and SEM of the tumor volume.
# Assemble the resulting series into a single summary DataFrame.
means = clean_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].mean()
medians = clean_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].median()
variances = clean_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].var()
stdevs = clean_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].std()
sems = clean_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].sem()
summary_table = pd.DataFrame({"Mean Tumor Volume":means,
                              "Median Tumor Volume":medians,
                              "Tumor Volume Variance":variances,
                              "Tumor Volume Std. Dev.":stdevs,
                              "Tumor Volume Std. Err.":sems})
summary_table


# ## Bar and Pie Charts

# In[175]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.
counts = clean_data['Drug Regimen'].value_counts()
counts.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("# of Observed mouse Timepoints")
plt.show()


# In[177]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.

import matplotlib.pyplot as plt

counts = clean_data['Drug Regimen'].value_counts()
counts.plot(kind="bar")

x = ['Drug_Regimen']
y = ['# of Observed mouse Timepoints']


plt.xlabel('Drug Regimen')
plt.xticks(rotation=90)
plt.ylabel('#of Observed mouse Timepoints')
plt.show()


# In[179]:


# Generate a pie chart, using Pandas, showing the distribution of unique female versus male mice used in the study
counts = clean_data.Sex.value_counts()
counts.plot(kind="pie",autopct='%.1f%%')
plt.show()


# In[181]:


# Generate a pie chart, using pyplot, showing the distribution of unique female versus male mice used in the study

import matplotlib.pyplot as plt

counts = clean_data.Sex.value_counts()
counts.plot(kind="pie",autopct='%.1f%%')
plt.show()



# Make the pie chart


# ## Quartiles, Outliers and Boxplots

# In[184]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:
# Capomulin, Ramicane, Infubinol, and Ceftamin

last_tumor = clean_data.groupby(["Mouse ID"])['Timepoint'].max()
last_tumor = last_tumor.reset_index()
print(last_tumor)


# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
merged_data_result = last_tumor.merge(clean_data, on=['Mouse ID','Timepoint'], how="left")
print(merged_data_result)


# In[186]:


# Put treatments into a list for for loop (and later for plot labels
treatment_list = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]


# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_list = []


# Calculate the IQR and quantitatively determine if there are any potential outliers.
for drug in treatment_list:

    # Locate the rows which contain mice on each drug and get the tumor volumes
    final_tumor_vol = merged_data_result.loc[merged_data_result["Drug Regimen"] == drug, 'Tumor Volume (mm3)']

    # add subset
    tumor_vol_list.append(final_tumor_vol)


    # Determine outliers using upper and lower bounds
    quartiles = final_tumor_vol.quantile([.25,.5,.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq-lowerq
    lower_bound = lowerq - (1.5*iqr)
    upper_bound = upperq + (1.5*iqr)
    outliers = final_tumor_vol.loc[(final_tumor_vol < lower_bound) | (final_tumor_vol > upper_bound)]
    print(f"{drug}'s potential outliers: {outliers}")


# In[188]:


# Generate a box plot that shows the distribution of the tumor volume for each treatment group.
mark_outlier_red = dict(markerfacecolor='red',markersize=12)
plt.boxplot(tumor_vol_list, labels = treatment_list, flierprops=mark_outlier_red)
plt.ylabel('Final Tumor Volum (mm3)')
plt.show()


# ## Line and Scatter Plots

# In[191]:


# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin
capomulin_table = clean_data.loc[clean_data['Drug Regimen'] == "Capomulin"]
mousedata = capomulin_table.loc[capomulin_table['Mouse ID']== '1509']

#create the line plot
plt.plot(mousedata['Timepoint'],mousedata['Tumor Volume (mm3)'])

#add labels and title
plt.xlabel('Timepoint (days)') 
plt.ylabel('Tumor Volume (mm3)') 
plt.title('Capomulin treatment of mouse 1509')

#Display the plot
plt.show()


# In[193]:


# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen
capomulin_table = clean_data.loc[clean_data['Drug Regimen'] == "Capomulin"]
capomulin_average = capomulin_table.groupby(['Mouse ID'])[['Weight (g)', 'Tumor Volume (mm3)']].mean()
plt.scatter(capomulin_average['Weight (g)'],capomulin_average['Tumor Volume (mm3)'])
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')
plt.show()


# ## Correlation and Regression

# In[196]:


# Calculate the correlation coefficient and a linear regression model
# for mouse weight and average observed tumor volume for the entire Capomulin regimen
corr=round(st.pearsonr(capomulin_average['Weight (g)'],capomulin_average['Tumor Volume (mm3)'])[0],2)
print(f"The correlation between mouse weight and the average tumor volume is {corr}")
model = st.linregress(capomulin_average['Weight (g)'],capomulin_average['Tumor Volume (mm3)'])

y_values = capomulin_average['Weight (g)']*model[0]+model[1]
plt.scatter(capomulin_average['Weight (g)'],capomulin_average['Tumor Volume (mm3)'])
plt.plot(capomulin_average['Weight (g)'],y_values,color="Orange")
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:




