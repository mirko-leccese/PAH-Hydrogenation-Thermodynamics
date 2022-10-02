# PAH-Hydrogenation-Thermodynamics

## Description
Fortran code to investigate the Thermodynamics of the stepwise-hydrogenation reaction of PAH molecules (such as coronene). TEST folder contains sample code outputs for coronene hydrogenation reaction, i.e. the stepwise addition of 24 H atoms. 

The code allows performing four types of analysis:
1. *Temperature analysis*: it computes partition functions for any hydrogenated species in a provided range of temperatures;
2. *Pressure analysis*: it computes pressure-dependent contribution to the Gibbs Free Energy in a provided range of partial pressures;
3. *Random Walk mode*: it simulates a mixture Free energy optimization starting from a promt composition of the mixture. The optimization details are given in the attached documentation;
4. *Phase diagram mode*: it simulates a phase diagram in a provided range of pressures ad temperatures.

## Requirements
Linux/MacOS operating system
gfortran compiler

## Usage
The code requires in the working directory all the files contaning the relevant thermodynamical data for each PAH specie. Such files have to be named as ```##-thermo.dat``` where "##" is the number of extra hydrogens attached to the PAH molecule (e.g. "00-thermo.dat" for the bare molecule, "01-thermo.dat" for the 1st hydrogenated specie, etc.). These files must have the following structure:

```
-921.283144 ! DFT Energy (eV)
78 78       ! alpha/beta electrons
300.09390   ! molar mass (amu)
1.          ! symmetry number
0.01606 0.01606 0.00803 ! Rotational temperatures (K)
! Harmonic modes (cm^-1)
87.5924 88.2717 121.7439 165.7634 226.1180 297.7808 298.5036 302.5380 302.7834 369.6797 369.8285 384.9954 385.3081 386.6606 457.4319 457.9992 478.4897 492.8873 494.8272 495.4836 533.2152 544.7634 545.1581 551.5028 565.5742 637.5588 668.5139 669.3492 679.9387 689.0734 689.2637 702.1300 741.6178 741.7717 781.7364 781.8213 782.6313 817.8555 819.7395 835.3841 835.5792 864.9778 865.5040 883.5007 931.4867 981.7058 989.5071 990.4130 1000.0023 1001.5280 1006.4376 1019.1244 1019.1841 1060.8544 1146.0668 1162.4495 1163.0700 1179.3256 1182.4624 1183.6040 1214.6503 1243.7384 1243.7411 1252.8409 1253.0290 1258.7816 1258.8983 1347.3134 1347.3513 1396.8842 1411.2007 1442.2085 1442.8080 1446.5947 1463.6857 1465.3024 1490.2839 1490.5021 1507.3945 1507.8301 1558.0811 1558.2463 1570.6701 1598.6724 1633.1635 1682.3940 1698.3566 1699.3847 1705.3554 1705.4376 3202.7136 3203.1322 3203.2432 3203.7320 3203.8366 3203.8697 3219.7280 3220.0861 3220.2414 3220.8791 3221.0400 3221.7570
```

