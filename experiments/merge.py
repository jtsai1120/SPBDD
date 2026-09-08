"""Merge the SPBDD probe log and the codeDistance tool log into one table."""
import re
import sys

spbdd_log, tools_log = sys.argv[1], sys.argv[2]

# SPBDD:  [[20,3,6]]     non-CSS       6          6       35.259    1785854
spbdd = {}
for line in open(spbdd_log):
    m = re.match(r"(\[\[\d+,\d+,\d+\]\])\s+(\S+)\s+(\d+)\s+(-?\d+)\s+([\d.]+)\s+(\d+)", line)
    if m:
        spbdd[m.group(1)] = {"type": m.group(2), "d": int(m.group(3)),
                             "got": int(m.group(4)), "secs": float(m.group(5)),
                             "nodes": int(m.group(6))}

# tools:  [[20,3,6]]     non-CSS     6 BZDistMW                 6      0.041
tools, methods = {}, []
for line in open(tools_log):
    m = re.match(r"(\[\[\d+,\d+,\d+\]\])\s+(\S+)\s+(\d+)\s+(\S+)\s+(-?\d+)\s+([\d.]+)\s*(.*)", line)
    if m:
        code, method, got, secs, note = m.group(1), m.group(4), int(m.group(5)), float(m.group(6)), m.group(7).strip()
        tools.setdefault(code, {})[method] = (got, secs, note)
        if method not in methods:
            methods.append(method)

def n_of(code):
    return int(code[2:].split(",")[0])

print(f"{'code':<13} {'type':<8} {'d':>3} | {'SPBDD':>10} | " +
      " | ".join(f"{m:>18}" for m in methods))
print("-" * (13 + 9 + 4 + 3 + 11 + 3 + 21 * len(methods)))

for code in sorted(set(spbdd) | set(tools), key=n_of):
    s = spbdd.get(code)
    row = f"{code:<13} "
    row += f"{(s['type'] if s else tools[code] and '?'):<8} " if s else f"{'':<8} "
    row += f"{(s['d'] if s else 0):>3} | "
    if s:
        flag = "" if s["got"] == s["d"] else "!"
        row += f"{s['secs']:>9.3f}{flag} | "
    else:
        row += f"{'--':>10} | "
    cells = []
    for m in methods:
        if code in tools and m in tools[code]:
            got, secs, note = tools[code][m]
            if "timeout" in note:
                cells.append(f"{'timeout':>18}")
            elif note:
                cells.append(f"{secs:>13.3f} {note[:4]:>4}")
            else:
                cells.append(f"{secs:>18.3f}")
        else:
            cells.append(f"{'--':>18}")
    print(row + " | ".join(cells))
