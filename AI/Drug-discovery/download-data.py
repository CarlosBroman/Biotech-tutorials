import pandas as pd
from chembl_webresource_client.new_client import new_client

# Create a new search
target = new_client.target
target_query = target.search('coronavirus')
targets = pd.DataFrame.from_dict(target_query)

# Select the desired target
selected_target = targets.target_chembl_id[6]
print(selected_target)

# Find the activity
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(res)
print(df.head(3))

df.standard_type.unique()

df.to_csv('bioactivity_data.csv', index = False)

df2 = df[df.standard_value.notna()]
print(df2)

bioactivity_class = []
for i in df2.standard_value:
    if float(i) >= 10000:
        bioactivity_class.append("inactive")
    elif float(i) <= 1000:
        bioactivity_class.append("active")
    else:
        bioactivity_class.append("intermediate")
            
selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
df3 = df2[selection]
df3 = pd.concat([df3, pd.DataFrame({'bioactivity_class':bioactivity_class})], axis=1)
print(df3)

df3.to_csv('bioactivity_preprocessed_data.csv', index=False)