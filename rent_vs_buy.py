import matplotlib.pyplot as plt

def calculate_rent_cost(monthly_rent, years, annual_increase_rate):
    total_rent_cost = 0
    for year in range(1, years + 1):
        total_rent_cost += monthly_rent * 12
        monthly_rent += monthly_rent * (annual_increase_rate / 100)
    return total_rent_cost

def calculate_buy_cost(home_price, down_payment_percent, loan_term_years, annual_interest_rate, property_tax_rate, maintenance_cost_percent):
    down_payment = home_price * (down_payment_percent / 100)
    loan_amount = home_price - down_payment
    monthly_interest_rate = annual_interest_rate / 12 / 100
    number_of_payments = loan_term_years * 12

    # Monthly mortgage payment calculation using the formula for an amortizing loan
    monthly_mortgage_payment = loan_amount * (monthly_interest_rate * (1 + monthly_interest_rate) ** number_of_payments) / ((1 + monthly_interest_rate) ** number_of_payments - 1)

    total_mortgage_cost = monthly_mortgage_payment * number_of_payments
    total_property_tax = home_price * (property_tax_rate / 100) * loan_term_years
    total_maintenance_cost = home_price * (maintenance_cost_percent / 100) * loan_term_years

    return down_payment + total_mortgage_cost + total_property_tax + total_maintenance_cost

def calculate_cumulative_buy_cost(home_price, down_payment_percent, loan_term_years, annual_interest_rate, property_tax_rate, maintenance_cost_percent, years):
    down_payment = home_price * (down_payment_percent / 100)
    loan_amount = home_price - down_payment
    monthly_interest_rate = annual_interest_rate / 12 / 100
    number_of_payments = loan_term_years * 12

    # Monthly mortgage payment calculation using the formula for an amortizing loan
    monthly_mortgage_payment = loan_amount * (monthly_interest_rate * (1 + monthly_interest_rate) ** number_of_payments) / ((1 + monthly_interest_rate) ** number_of_payments - 1)

    cumulative_costs = []
    total_mortgage_cost = 0
    for year in range(1, years + 1):
        if year <= loan_term_years:
            total_mortgage_cost += monthly_mortgage_payment * 12
        total_property_tax = home_price * (property_tax_rate / 100) * year
        total_maintenance_cost = home_price * (maintenance_cost_percent / 100) * year
        cumulative_cost = down_payment + total_mortgage_cost + total_property_tax + total_maintenance_cost
        cumulative_costs.append(cumulative_cost)

    return cumulative_costs

def main():
    # Inputs
    monthly_rent = 1500
    annual_rent_increase_rate = 2  # Annual increase rate for rent
    home_price = 300000
    down_payment_percent = 20
    loan_term_years = 20
    annual_interest_rate = 3.5
    property_tax_rate = 1.2
    maintenance_cost_percent = 1

    # Calculations
    years = list(range(1, 31))
    rent_costs = [calculate_rent_cost(monthly_rent, year, annual_rent_increase_rate) for year in years]
    cumulative_buy_costs = calculate_cumulative_buy_cost(home_price, down_payment_percent, loan_term_years, annual_interest_rate, property_tax_rate, maintenance_cost_percent, max(years))
    buy_costs = cumulative_buy_costs[:len(years)]

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(years, rent_costs, label='Rent Cost')
    plt.plot(years, cumulative_buy_costs, label='Cumulative Buy Cost')
    plt.xlabel('Years')
    plt.ylabel('Total Cost ($)')
    plt.title('Rent vs Cumulative Buy Cost Comparison')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()