import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import sympy as sp

class SeriesRLCircuit:
    def __init__(self, resistance, inductance):
        """
        Initialize the series RL circuit with resistance and inductance values.
        
        Args:
            resistance (float): Resistance value in ohms
            inductance (float): Inductance value in henries
        """
        self.R = resistance
        self.L = inductance
        
    def differential_equation(self, t, y, voltage_func):
        """
        Define the differential equation for the RL circuit: L*di/dt + R*i = v(t)
        Rearranged to: di/dt = (v(t) - R*i) / L
        
        Args:
            t (float): Time in seconds
            y (list): Current state [i]
            voltage_func (function): Function that returns voltage at time t
            
        Returns:
            list: Derivative of the state [di/dt]
        """
        i = y[0]
        v = voltage_func(t)
        di_dt = (v - self.R * i) / self.L
        return [di_dt]
    
    def simulate(self, voltage_func, t_span, num_points=1000):
        """
        Simulate the circuit response over time.
        
        Args:
            voltage_func (function): Function that returns voltage at time t
            t_span (tuple): Time span (t_start, t_end) in seconds
            num_points (int): Number of points to evaluate
            
        Returns:
            tuple: (time_points, currents, resistor_voltages)
        """
        # Initial condition: zero current
        y0 = [0]
        
        # Solve the differential equation
        t_eval = np.linspace(t_span[0], t_span[1], num_points)
        solution = solve_ivp(
            fun=lambda t, y: self.differential_equation(t, y, voltage_func),
            t_span=t_span,
            y0=y0,
            t_eval=t_eval,
            method='RK45'
        )
        
        # Extract the current values
        currents = solution.y[0]
        
        # Calculate voltage across the resistor: V_R = I * R
        resistor_voltages = currents * self.R
        
        return solution.t, currents, resistor_voltages

# Predefined waveform functions
def sine_wave(amplitude, frequency, phase=0):
    """Return a function that generates a sine wave."""
    return lambda t: amplitude * np.sin(2 * np.pi * frequency * t + phase)

def square_wave(amplitude, frequency, duty_cycle=0.5):
    """Return a function that generates a square wave."""
    return lambda t: amplitude if (t * frequency) % 1 < duty_cycle else 0

def sawtooth_wave(amplitude, frequency):
    """Return a function that generates a sawtooth wave."""
    return lambda t: amplitude * ((t * frequency) % 1)

def triangle_wave(amplitude, frequency):
    """Return a function that generates a triangle wave."""
    return lambda t: amplitude * (2 * abs(2 * ((t * frequency) % 1) - 1) - 1)

def dc_voltage(amplitude):
    """Return a function that generates a constant DC voltage."""
    return lambda t: amplitude

def custom_waveform(expression_str, t_symbol='t'):
    """
    Create a custom waveform from a mathematical expression string.
    
    Args:
        expression_str (str): Mathematical expression with 't' as the time variable
        t_symbol (str): Symbol to use for time variable
        
    Returns:
        function: A function that evaluates the expression at time t
    """
    # Create a symbolic variable for time
    t = sp.Symbol(t_symbol)
    
    # Parse the expression
    expr = sp.sympify(expression_str)
    
    # Create a lambda function
    func = sp.lambdify(t, expr, modules=["numpy"])
    
    return func

def main():
    print("\n===== Series RL Circuit Simulator =====\n")
    
    # Get circuit parameters
    while True:
        try:
            resistance = float(input("Enter resistance (R) in ohms: "))
            if resistance <= 0:
                raise ValueError("Resistance must be positive")
            break
        except ValueError as e:
            print(f"Invalid input: {e}")
    
    while True:
        try:
            inductance = float(input("Enter inductance (L) in henries: "))
            if inductance <= 0:
                raise ValueError("Inductance must be positive")
            break
        except ValueError as e:
            print(f"Invalid input: {e}")
    
    # Create the circuit
    circuit = SeriesRLCircuit(resistance, inductance)
    
    # Get waveform type
    print("\nSelect voltage source waveform:")
    print("1. Sine wave")
    print("2. Square wave")
    print("3. Sawtooth wave")
    print("4. Triangle wave")
    print("5. DC voltage")
    print("6. Custom waveform")
    
    while True:
        try:
            choice = int(input("\nEnter your choice (1-6): "))
            if choice not in range(1, 7):
                raise ValueError("Please enter a number between 1 and 6")
            break
        except ValueError as e:
            print(f"Invalid input: {e}")
    
    # Get waveform parameters
    if choice in [1, 2, 3, 4]:
        while True:
            try:
                amplitude = float(input("Enter voltage amplitude (V): "))
                frequency = float(input("Enter frequency (Hz): "))
                if frequency <= 0:
                    raise ValueError("Frequency must be positive")
                break
            except ValueError as e:
                print(f"Invalid input: {e}")
        
        if choice == 1:  # Sine wave
            voltage_func = sine_wave(amplitude, frequency)
            waveform_name = f"Sine wave ({amplitude}V, {frequency}Hz)"
        elif choice == 2:  # Square wave
            duty_cycle = float(input("Enter duty cycle (0-1): "))
            voltage_func = square_wave(amplitude, frequency, duty_cycle)
            waveform_name = f"Square wave ({amplitude}V, {frequency}Hz, duty cycle={duty_cycle})"
        elif choice == 3:  # Sawtooth wave
            voltage_func = sawtooth_wave(amplitude, frequency)
            waveform_name = f"Sawtooth wave ({amplitude}V, {frequency}Hz)"
        else:  # Triangle wave
            voltage_func = triangle_wave(amplitude, frequency)
            waveform_name = f"Triangle wave ({amplitude}V, {frequency}Hz)"
        
    elif choice == 5:  # DC voltage
        amplitude = float(input("Enter DC voltage (V): "))
        voltage_func = dc_voltage(amplitude)
        waveform_name = f"DC voltage ({amplitude}V)"
        frequency = 0  # No frequency for DC
        
    else:  # Custom waveform
        print("\nEnter a custom mathematical expression using 't' as the time variable.")
        print("Examples: '5*sin(2*pi*10*t)', 't**2', 'exp(-t)*sin(2*pi*5*t)'")
        print("Available functions: sin, cos, tan, exp, log, sqrt, abs, etc.")
        expression_str = input("\nExpression: ")
        
        try:
            # Test the expression
            custom_func = custom_waveform(expression_str)
            custom_func(0)  # Test with t=0
            
            voltage_func = custom_func
            waveform_name = f"Custom waveform: {expression_str}"
            
            frequency = float(input("Enter characteristic frequency for plotting (Hz): "))
        except Exception as e:
            print(f"Error in expression: {e}")
            print("Using default sine wave instead.")
            amplitude = 5.0
            frequency = 1.0
            voltage_func = sine_wave(amplitude, frequency)
            waveform_name = f"Sine wave ({amplitude}V, {frequency}Hz)"
    
    # Set simulation parameters
    if frequency > 0:
        # For oscillating signals, simulate a few periods
        periods_to_show = 5
        t_end = periods_to_show / frequency
    else:
        # For DC signals, show the transient response
        time_constant = circuit.L / circuit.R
        t_end = 5 * time_constant  # 5 time constants to reach steady state
    
    t_span = (0, t_end)
    
    # Run simulation
    print("\nRunning simulation...")
    times, currents, resistor_voltages = circuit.simulate(voltage_func, t_span)
                                                                            
    # Calculate and display results
    max_voltage = np.max(resistor_voltages)
    min_voltage = np.min(resistor_voltages)
    rms_voltage = np.sqrt(np.mean(resistor_voltages**2))
    
    print("\n===== Results =====")
    print(f"Time constant (L/R): {circuit.L/circuit.R:.6f} seconds")
    print(f"Maximum voltage across resistor: {max_voltage:.4f} V")
    print(f"Minimum voltage across resistor: {min_voltage:.4f} V")
    print(f"RMS voltage across resistor: {rms_voltage:.4f} V")
    
    
    # Plot the results
    while True:
        doIprint = input("Do you want to plot the results? y/n?").lower()
    
        if doIprint == "y":
            plt.figure(figsize=(12, 8))
            
            # Input voltage
            plt.subplot(3, 1, 1)
            input_voltages = [voltage_func(t) for t in times]
            plt.plot(times, input_voltages)
            plt.title(f"Input Voltage: {waveform_name}")
            plt.ylabel("Voltage (V)")
            plt.grid(True)
            
            # Current
            plt.subplot(3, 1, 2)
            plt.plot(times, currents)
            plt.title("Current through the Circuit")
            plt.ylabel("Current (A)")
            plt.grid(True)
            
            # Resistor voltage
            plt.subplot(3, 1, 3)
            plt.plot(times, resistor_voltages)
            plt.title("Voltage across Resistor")
            plt.xlabel("Time (s)")
            plt.ylabel("Voltage (V)")
            plt.grid(True)
            
            plt.tight_layout()
            plt.show()
            break
        elif doIprint == "n":
            print("Exit")
            break
        else:
            print("Invalid input!!. type 'y' or 'n'")
            continue
        

if __name__ == "__main__":
    main()