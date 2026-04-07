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

<img width="640" height="480" alt="plot-bias-R1-k31-abs" src="https://github.com/user-attachments/assets/a208e124-ddec-42cb-8419-3c4181b54043" />
<img width="640" height="480" alt="plot-bias-R1-k31-rel" src="https://github.com/user-attachments/assets/23c3e776-9e3d-4ab9-a389-5c4ce2c12380" />
