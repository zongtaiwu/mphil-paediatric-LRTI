import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.patches import Patch
# Read the CSV file
f1_scores = pd.read_csv('/Users/zongt1/VAP/project-main/results/Pred-host/F1_scores.csv')

# Continue with the rest of the code
melted_f1_scores = f1_scores.melt(var_name='Model', value_name='F1_Score')
palette = {
    '7-LR': 'firebrick', '14-LR': 'darkcyan', 'all-LR': 'gray',
    '7-RF': 'firebrick', '14-RF': 'darkcyan', 'all-RF': 'gray',
    '7-SVM': 'firebrick', '14-SVM': 'darkcyan', 'all-SVM': 'gray',
    '7-KNN': 'firebrick', '14-KNN': 'darkcyan', 'all-KNN': 'gray',
    '7-XGBoost': 'firebrick', '14-XGBoost': 'darkcyan', 'all-XGBoost': 'gray'
}

order = ['7-LR', '14-LR', 'all-LR',
         '7-RF', '14-RF', 'all-RF',
         '7-SVM', '14-SVM', 'all-SVM',
         '7-KNN', '14-KNN', 'all-KNN',
         '7-XGBoost', '14-XGBoost', 'all-XGBoost']

plt.figure(figsize=(18, 10))  # Figure size
sns.boxplot(data=melted_f1_scores, x='Model', y='F1_Score', palette=palette, order=order)
plt.xticks(rotation=30, fontsize=18)  # Increase the font size of x-axis labels
plt.yticks(fontsize=18)  # Increase the font size of y-axis labels
plt.title('F1 Scores for Different Models', fontsize=20, pad = 80)  # Increase the title font size

# Add vertical lines to separate the groups
plt.axvline(x=2.5, color='black', linestyle='--')
plt.axvline(x=5.5, color='black', linestyle='--')
plt.axvline(x=8.5, color='black', linestyle='--')
plt.axvline(x=11.5, color='black', linestyle='--')

# Remove the upper boundary (spine)
plt.gca().spines['top'].set_visible(False)


legend_elements = [
    Patch(facecolor='firebrick', edgecolor='black', label='7 genes'),
    Patch(facecolor='darkcyan', edgecolor='black', label='14 genes'),
    Patch(facecolor='gray', edgecolor='black', label='All genes')
]

# Add the legend to the plot with increased font size
plt.legend(handles=legend_elements, loc='lower right', fontsize=18)

plt.show()


plt.savefig('/Users/zongt1/VAP/project-main/results/Pred/boxplot.png')