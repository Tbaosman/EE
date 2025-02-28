import numpy as np
import matplotlib.pyplot as plt
from math import pi

def calculate_resistor_voltage(supply_voltage, resistance, inductance, frequency):
    """
    Calculate the voltage across a resistor in a series RL circuit.
    
    Parameters:
    supply_voltage (float): The voltage of the power supply in Volts
    resistance (float): The resistance in Ohms
    inductance (float): The inductance in Henries
    frequency (float): The frequency in Hertz
    
    Returns:
    float: The voltage across the resistor in Volts
    """
    # Calculate angular frequency
    omega = 2 * pi * frequency
    
    # Calculate the impedance of the inductor
    inductive_reactance = omega * inductance
    
    # Calculate the total impedance of the circuit
    total_impedance = np.sqrt(resistance**2 + inductive_reactance**2)
    
    # Calculate the phase angle
    phase_angle = np.arctan(inductive_reactance / resistance)
    
    # Calculate the current through the circuit
    current = supply_voltage / total_impedance
    
    # Calculate the voltage across the resistor
    resistor_voltage = current * resistance
    
    return resistor_voltage, phase_angle

def main():
    print("Series RL Circuit Calculator")
    print("===========================")
    
    # Get user inputs
    supply_voltage = float(input("Enter the supply voltage (V): "))
    resistance = float(input("Enter the resistance (Ω): "))
    inductance = float(input("Enter the inductance (H): "))
    frequency = float(input("Enter the frequency (Hz): "))
    
    # Calculate the voltage across the resistor
    resistor_voltage, phase_angle = calculate_resistor_voltage(supply_voltage, resistance, inductance, frequency)
    
    # Display the results
    print("\nResults:")
    print(f"Voltage across the resistor: {resistor_voltage:.2f} V")
    print(f"Current phase angle: {phase_angle * 180/pi:.2f} degrees")
    
    # Calculate impedance
    omega = 2 * pi * frequency
    inductive_reactance = omega * inductance
    total_impedance = np.sqrt(resistance**2 + inductive_reactance**2)
    print(f"Total circuit impedance: {total_impedance:.2f} Ω")
    print(f"Current through the circuit: {supply_voltage/total_impedance:.4f} A")
    
# Visualization options
    print("\nVisualization Options:")
    print("1. Resistor voltage vs. frequency")
    print("2. Phase angle vs. frequency")
    print("3. Impedance vs. frequency")
    print("4. Voltage division (resistor/inductor) vs. frequency")
    print("5. Time domain voltage and current waveforms")
    print("6. Phasor diagram at current frequency")
    print("7. Bode plot (magnitude and phase)")
    print("8. Exit")
    
    while True:
        choice = input("\nSelect a visualization (1-8): ")
        
        if choice == '1':
            plot_frequency_response(supply_voltage, resistance, inductance, "resistor")
        elif choice == '2':
            plot_phase_response(resistance, inductance)
        elif choice == '3':
            plot_impedance_response(resistance, inductance)
        elif choice == '4':
            plot_voltage_division(supply_voltage, resistance, inductance)
        elif choice == '5':
            plot_time_domain(supply_voltage, resistance, inductance, frequency)
        elif choice == '6':
            plot_phasor_diagram(supply_voltage, resistance, inductance, frequency)
        elif choice == '7':
            plot_bode(resistance, inductance)
        elif choice == '8':
            break
        else:
            print("Invalid option. Please try again.")

def plot_frequency_response(supply_voltage, resistance, inductance, component_type):
    frequencies = np.logspace(0, 5, 1000)  # from 1 Hz to 100 kHz
    resistor_voltages = []
    inductor_voltages = []
    
    for freq in frequencies:
        r_v, l_v, _, _ = calculate_resistor_voltage(supply_voltage, resistance, inductance, freq)
        resistor_voltages.append(r_v)
        inductor_voltages.append(l_v)
    
    plt.figure(figsize=(10, 6))
    if component_type == "resistor":
        plt.semilogx(frequencies, resistor_voltages)
        plt.ylabel('Resistor Voltage (V)')
        plt.title('Voltage across Resistor vs. Frequency')
    else:
        plt.semilogx(frequencies, inductor_voltages)
        plt.ylabel('Inductor Voltage (V)')
        plt.title('Voltage across Inductor vs. Frequency')
    
    plt.xlabel('Frequency (Hz)')
    plt.grid(True)
    plt.show()

def plot_phase_response(resistance, inductance):
    frequencies = np.logspace(0, 5, 1000)  # from 1 Hz to 100 kHz
    phases = []
    
    for freq in frequencies:
        omega = 2 * pi * freq
        inductive_reactance = omega * inductance
        phase = np.arctan(inductive_reactance / resistance)
        phases.append(phase * 180 / pi)  # Convert to degrees
    
    plt.figure(figsize=(10, 6))
    plt.semilogx(frequencies, phases)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase Angle (degrees)')
    plt.title('Phase Angle vs. Frequency')
    plt.grid(True)
    plt.show()

def plot_impedance_response(resistance, inductance):
    frequencies = np.logspace(0, 5, 1000)  # from 1 Hz to 100 kHz
    impedances = []
    r_components = []
    l_components = []
    
    for freq in frequencies:
        omega = 2 * pi * freq
        inductive_reactance = omega * inductance
        impedance = np.sqrt(resistance**2 + inductive_reactance**2)
        impedances.append(impedance)
        r_components.append(resistance)
        l_components.append(inductive_reactance)
    
    plt.figure(figsize=(10, 6))
    plt.semilogx(frequencies, impedances, label='Total Impedance')
    plt.semilogx(frequencies, r_components, label='Resistance', linestyle='--')
    plt.semilogx(frequencies, l_components, label='Inductive Reactance', linestyle=':')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Impedance (Ω)')
    plt.title('Impedance Components vs. Frequency')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_voltage_division(supply_voltage, resistance, inductance):
    frequencies = np.logspace(0, 5, 1000)  # from 1 Hz to 100 kHz
    vr_ratio = []
    vl_ratio = []
    
    for freq in frequencies:
        r_v, l_v, _, _ = calculate_resistor_voltage(supply_voltage, resistance, inductance, freq)
        vr_ratio.append(r_v / supply_voltage)
        vl_ratio.append(l_v / supply_voltage)
    
    plt.figure(figsize=(10, 6))
    plt.semilogx(frequencies, vr_ratio, label='Resistor Voltage Ratio')
    plt.semilogx(frequencies, vl_ratio, label='Inductor Voltage Ratio')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Voltage Ratio (V/Vsupply)')
    plt.title('Voltage Division vs. Frequency')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_time_domain(supply_voltage, resistance, inductance, frequency):
    omega = 2 * pi * frequency
    inductive_reactance = omega * inductance
    impedance = np.sqrt(resistance**2 + inductive_reactance**2)
    phase_angle = np.arctan(inductive_reactance / resistance)
    current_amplitude = supply_voltage / impedance
    
    # Create time points for one cycle
    t = np.linspace(0, 1/frequency, 1000)
    
    # Calculate voltages and currents
    supply_v = supply_voltage * np.sin(omega * t)
    current = current_amplitude * np.sin(omega * t - phase_angle)
    resistor_v = resistance * current
    inductor_v = supply_v - resistor_v
    
    # Plot
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 1, 1)
    plt.plot(t, supply_v, label='Supply Voltage')
    plt.plot(t, resistor_v, label='Resistor Voltage')
    plt.plot(t, inductor_v, label='Inductor Voltage')
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    plt.title('Voltage Waveforms')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(2, 1, 2)
    plt.plot(t, current, label='Current')
    plt.xlabel('Time (s)')
    plt.ylabel('Current (A)')
    plt.title('Current Waveform')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()

def plot_phasor_diagram(supply_voltage, resistance, inductance, frequency):
    omega = 2 * pi * frequency
    inductive_reactance = omega * inductance
    impedance = np.sqrt(resistance**2 + inductive_reactance**2)
    phase_angle = np.arctan(inductive_reactance / resistance)
    current_amplitude = supply_voltage / impedance
    
    resistor_voltage = current_amplitude * resistance
    inductor_voltage = current_amplitude * inductive_reactance
    
    # Create figure
    plt.figure(figsize=(8, 8))
    ax = plt.subplot(111, projection='polar')
    
    # Plot phasors
    # Current at 0 degrees (reference)
    ax.arrow(0, 0, 0, current_amplitude, alpha=0.5, width=0.02, 
            edgecolor='blue', facecolor='blue', lw=2, zorder=5)
    
    # Supply voltage
    ax.arrow(0, 0, phase_angle, supply_voltage, alpha=0.5, width=0.02,
            edgecolor='black', facecolor='black', lw=2, zorder=5)
    
    # Resistor voltage (in phase with current)
    ax.arrow(0, 0, 0, resistor_voltage, alpha=0.5, width=0.02,
            edgecolor='red', facecolor='red', lw=2, zorder=5)
    
    # Inductor voltage (90 degrees ahead of current)
    ax.arrow(0, 0, pi/2, inductor_voltage, alpha=0.5, width=0.02,
            edgecolor='green', facecolor='green', lw=2, zorder=5)
    
    # Add labels
    ax.text(0, current_amplitude/2, "Current", color='blue', ha='right', va='bottom')
    ax.text(phase_angle, supply_voltage/2, "Supply V", color='black', ha='left', va='bottom')
    ax.text(0, resistor_voltage/2, "Resistor V", color='red', ha='left', va='bottom')
    ax.text(pi/2, inductor_voltage/2, "Inductor V", color='green', ha='left', va='bottom')
    
    ax.set_title('Phasor Diagram for RL Circuit')
    plt.show()

def plot_bode(resistance, inductance):
    frequencies = np.logspace(0, 5, 1000)  # from 1 Hz to 100 kHz
    magnitude_r = []
    magnitude_l = []
    phases = []
    
    for freq in frequencies:
        omega = 2 * pi * freq
        inductive_reactance = omega * inductance
        impedance = np.sqrt(resistance**2 + inductive_reactance**2)
        phase = np.arctan(inductive_reactance / resistance)
        
        # Magnitude in dB
        magnitude_r.append(20 * np.log10(resistance / impedance))
        magnitude_l.append(20 * np.log10(inductive_reactance / impedance))
        phases.append(phase * 180 / pi)  # Convert to degrees
    
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 1, 1)
    plt.semilogx(frequencies, magnitude_r, label='Resistor')
    plt.semilogx(frequencies, magnitude_l, label='Inductor')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (dB)')
    plt.title('Bode Plot - Magnitude')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(2, 1, 2)
    plt.semilogx(frequencies, phases)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase (degrees)')
    plt.title('Bode Plot - Phase')
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
