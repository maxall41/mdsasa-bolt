def save(analysis, filename):
    with open("bench/" + filename, "w") as f:
        for frame in analysis.results.total_area:
            # frame is float
            f.write(f"{frame}\n")
