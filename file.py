import pandas as pd

# Read the large CSV file
large_csv_file = "temp_route_data.csv"

# Read the first 100 rows
df_first_100_rows = pd.read_csv(large_csv_file, nrows=10)

# Save the first 100 rows to a new CSV file
small_csv_file = "temp_route_data_small).csv"
df_first_100_rows.to_csv(small_csv_file, index=False)

print(f"The first 10 rows have been written to {small_csv_file}")
