Intelligent Reflecting Surface Operation under Predictable Receiver Mobility: A Continuous Time Propagation Model
==================

This code package is related to the following scientific article:

Bho Matthiesen, Emil Björnson, Elisabeth De Carvalho, and Petar Popovski, "[Intelligent Reflecting Surface Operation under Predictable Receiver Mobility: A Continuous Time Propagation Model](https://doi.org/10.1109/LWC.2020.3024781)," IEEE Wireless Communications Letters.


## Abstract of Article

The operation of an intelligent reflecting surface (IRS) under predictable receiver mobility is investigated. We develop a continuous time system model for multipath channels and discuss the optimal IRS configuration with respect to received power, Doppler spread, and delay spread. It is shown that the received power can be maximized without adding Doppler spread to the system. In a numerical case study, we show that an IRS having the size of just two large billboards can improve the link budget of ground to Low Earth Orbit (LEO) satellite links by up to 6 dB. It also adds a second, almost equivalently strong, communication path that improves the link reliability.


## Requirements & Usage

This code was tested with [Python](https://python.org) 3.8.5, [NumPy](https://numpy.org) 1.19.1, [SciPy](https://scipy.org) 1.5.2, [Matplotlib](https://matplotlib.org) 3.3.1, and [pandas](https://pandas.pydata.org) 1.1.1.

Run the script as `python3 satcom.py` to reproduce Figures 3 and 4 in the article cited above. Consider changing the parameters `MC` and `Tpoints` in lines 315-7 to reduce the computation time.

## Acknowledgements

This research as supported in part by the German Research Foundation (DFG) under Germany's Excellence Strategy (EXC 2077 at University of Bremen, University Allowance), in part by ELLIIT and in part by the Danish Council for Independent Research DFF-701700271.


## License and Referencing

This program is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

