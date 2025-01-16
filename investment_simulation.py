import matplotlib.pyplot as plt

# Define the investment parameters
initial_investment = 10000  # Initial investment amount in dollars
annual_interest_rate = 0.05  # Annual interest rate (5%)
annual_inflation_rate = 0.02  # Annual inflation rate (2%)
annual_interest_rate_1 = 0.02  # Annual interest rate (2%)
annual_inflation_rate_1 = 0.02  # Annual inflation rate (2%)
annual_interest_rate_2 = 0.00  # Annual interest rate (0%)
annual_inflation_rate_2 = 0.02  # Annual inflation rate (2%)
years = 30  # Number of years

# Calculate the investment value for each year
investment_values = [initial_investment]
inflation_adjusted_values = [initial_investment]
investment_values_1 = [initial_investment]
inflation_adjusted_values_1 = [initial_investment]
investment_values_2 = [initial_investment]
inflation_adjusted_values_2 = [initial_investment]
for year in range(1, years + 1):
    investment_value = investment_values[-1] * (1 + annual_interest_rate)
    inflation_adjusted_value = inflation_adjusted_values[-1] * (1 + annual_interest_rate - annual_inflation_rate)
    investment_values.append(investment_value)
    inflation_adjusted_values.append(inflation_adjusted_value)

    investment_value_1 = investment_values_1[-1] * (1 + annual_interest_rate_1)
    inflation_adjusted_value_1 = inflation_adjusted_values_1[-1] * (1 + annual_interest_rate_1 - annual_inflation_rate_1)
    investment_values_1.append(investment_value_1)
    inflation_adjusted_values_1.append(inflation_adjusted_value_1)

    investment_value_2 = investment_values_2[-1] * (1 + annual_interest_rate_2)
    inflation_adjusted_value_2 = inflation_adjusted_values_2[-1] * (1 + annual_interest_rate_2 - annual_inflation_rate_2)
    investment_values_2.append(investment_value_2)
    inflation_adjusted_values_2.append(inflation_adjusted_value_2)

# Generate the years list
years_list = list(range(0, years + 1))

# Plot the investment growth over time
plt.figure(figsize=(10, 6))
plt.plot(years_list, investment_values, marker='o', linestyle='-', color='b', label='Investment Value')
plt.plot(years_list, inflation_adjusted_values, marker='o', linestyle='--', color='r', label='Inflation Adjusted Value')
plt.plot(years_list, investment_values_1, marker='o', linestyle='-', color='g', label='Investment Value 1')
plt.plot(years_list, inflation_adjusted_values_1, marker='o', linestyle='--', color='y', label='Inflation Adjusted Value 1')
plt.plot(years_list, investment_values_2, marker='o', linestyle='-', color='m', label='Investment Value 2')
plt.plot(years_list, inflation_adjusted_values_2, marker='o', linestyle='--', color='c', label='Inflation Adjusted Value 2')
plt.title('Investment Growth Over Time')
plt.xlabel('Years')
plt.ylabel('Investment Value ($)')
plt.legend()
plt.grid(True)
plt.show()
    



