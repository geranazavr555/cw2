C:\Windows\system32\wsl.exe --distribution Ubuntu-22.04 --exec /bin/bash -c "export PARLAY_NUM_THREADS=4 && export CILK_NWORKERS=4 && cd /mnt/c/Programing/sem7/par-algo/cw2/cmake-build-release && /mnt/c/Programing/sem7/par-algo/cw2/cmake-build-release/cw2"

Current size: 50
sequential micros: 46038
  parallel micros: 23255
Boost: 1.97972
Current size: 100
sequential micros: 336043
  parallel micros: 155300
Boost: 2.16383
Current size: 250
sequential micros: 6529592
  parallel micros: 2360150
Boost: 2.7666
Current size: 500
sequential micros: 56995465
  parallel micros: 20329010
Boost: 2.80365

C:\Windows\system32\wsl.exe --distribution Ubuntu-22.04 --exec /bin/bash -c "export ASAN_SYMBOLIZER_PATH=/opt/opencilk/bin/llvm-symbolizer && export CILK_NWORKERS=8 && cd /mnt/c/Programing/sem7/par-algo/cw2/cmake-build-release && /mnt/c/Programing/sem7/par-algo/cw2/cmake-build-release/cw2"
Current size: 50
sequential micros: 62691
  parallel micros: 18219
Boost: 3.44097
Current size: 100
sequential micros: 480799
  parallel micros: 130982
Boost: 3.67071
Current size: 250
sequential micros: 8937556
  parallel micros: 1935507
Boost: 4.61768
Current size: 500
sequential micros: 77602337
  parallel micros: 16095265
Boost: 4.82144