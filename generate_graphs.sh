cargo flamegraph -o cycles.svg -c "record -e cycles -c 100 --call-graph dwarf -g" -- list.txt
cargo flamegraph -o instructions.svg -c "record -e instructions -c 100 --call-graph dwarf -g" -- list.txt
cargo flamegraph --flamechart -o cycles-chart.svg -c "record -e cycles -c 100 --call-graph dwarf -g" -- list.txt
cargo flamegraph --flamechart -o instructions-chart.svg -c "record -e instructions -c 100 --call-graph dwarf -g" -- list.txt

