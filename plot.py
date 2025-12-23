import csv
import matplotlib.pyplot as plt

def read_csv(path):
    data = {
        "n": [],
        "time": [],
        "muls": [],
        "adds": [],
    }
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data["n"].append(int(row["n"]))
            data["time"].append(float(row["mean_s"]))
            data["muls"].append(int(row["muls"]))
            data["adds"].append(int(row["adds"]))
    return data

naive = read_csv("naive.csv")
strassen = read_csv("strassen.csv")
alpha = read_csv("alphaevolve.csv")

plt.figure()
plt.plot(naive["n"], naive["time"], marker="o", label="naive")
plt.plot(strassen["n"], strassen["time"], marker="o", label="strassen")
plt.plot(alpha["n"], alpha["time"], marker="o", label="alphaevolve")
plt.xscale("log", base=2)
plt.xlabel("Matrix size (n)")
plt.ylabel("Mean time (seconds)")
plt.title("Mean execution time")
plt.legend()
plt.grid(True)
plt.savefig("time.jpg")
plt.close()

plt.figure()
plt.plot(naive["n"], naive["muls"], marker="o", label="naive")
plt.plot(strassen["n"], strassen["muls"], marker="o", label="strassen")
plt.plot(alpha["n"], alpha["muls"], marker="o", label="alphaevolve")
plt.xscale("log", base=2)
plt.xlabel("Matrix size (n)")
plt.ylabel("Number of multiplications")
plt.title("Multiplications count")
plt.legend()
plt.grid(True)
plt.savefig("muls.jpg")
plt.close()

plt.figure()
plt.plot(naive["n"], naive["adds"], marker="o", label="naive")
plt.plot(strassen["n"], strassen["adds"], marker="o", label="strassen")
plt.plot(alpha["n"], alpha["adds"], marker="o", label="alphaevolve")
plt.xscale("log", base=2)
plt.xlabel("Matrix size (n)")
plt.ylabel("Number of additions")
plt.title("Additions count")
plt.legend()
plt.grid(True)
plt.savefig("adds.jpg")
plt.close()
