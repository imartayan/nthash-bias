# NtHash bias

Visualize the bias of NtHash for consecutive k-mers, and how changing the rotation factor impacts it

```sh
# plot bias for R=1 and R=7 using 4 threads on a random seq
cargo r -r -- -r 1 7 -t 4 | python plot-bias.py

# passing -i <input.fa> will use a specific sequence instead
cargo r -r -- -r 1 7 -t 4 -i <input.fa> | python plot-bias.py

# output multiple plot formats (png by default)
cargo r -r -- -r 1 7 -t 4 | python plot-bias.py -f pdf svg
```
