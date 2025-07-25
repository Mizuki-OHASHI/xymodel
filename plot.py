from sys import argv
import matplotlib.pyplot as plt


def main():
    if len(argv) < 2:
        print("Usage: python plot.py <data_file>")
        return

    data_file = argv[1]
    with open(data_file, "r") as f:
        lines = f.readlines()

    data = {}
    label = None
    for line in lines:
        line = line.strip()
        if not line:
            continue
        elif line.startswith("<"):
            label = line[1:]
            data[label] = []
        elif line.endswith(">"):
            assert label == line[:-1]
            label = None
        elif "=" in line:
            keyvalue = line.split(" ")
            data[label].append({k: v for k, v in (kv.split("=") for kv in keyvalue)})
        else:
            data[label].append(line)

    betas = [float(item["beta"]) for item in data["helicity_modulus_results"]]
    helicity_moduli = [
        float(item["helicity_modulus"]) for item in data["helicity_modulus_results"]
    ]
    input_params = data["input_parameters"][0]
    input_params_str = ""
    MAX_LINE_LENGTH = 64
    current_length = 0
    for k, v in input_params.items():
        if k in ["N_SITES", "N_REPLICAS"]:
            continue
        if current_length + len(k) + len(v) + 3 > MAX_LINE_LENGTH:
            input_params_str = input_params_str.rstrip(", ") + "\n"
            current_length = 0
        kv = f"{k}={v}, "
        input_params_str += kv
        current_length += len(kv)
    input_params_str = input_params_str.rstrip(", ")

    plt.figure(figsize=(8, 6))
    plt.plot(betas, helicity_moduli, "bo")
    plt.title(
        f"Helicity Modulus of XY Model vs. Inverse Temperature (took {data['execution_time'][0]} seconds)\n{input_params_str}"
    )
    plt.xlabel("Inverse Temperature (Beta)")
    plt.ylabel("Helicity Modulus")
    plt.tight_layout()
    plt.savefig(f"figs/{data_file}.png", dpi=300)


if __name__ == "__main__":
    main()
