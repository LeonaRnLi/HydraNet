# From Interference Mitigation to Toleration: Pathway to Practical Spatial Reuse in LPWANs
HydraNet can improve spatial reuse and optimize spectrum efficiency in LPWANs.

This webpage contains source code to deploy HydraNet.


## Overview
Low Power Wide-Area Networks (LPWANs) suffer from sever interference due to their long range transmissions. Traditional spatial reuse solutions share a common philosophy, i.e., interference mitigation. 

HydraNet is a new paradigm to facilitate simultaneous LPWAN transmissions on the same frequency,fundamentally evolving the current philosophy of interference management—from mitigating interference to tolerating its presence.

The PDF file is in the main folder.

***Reference:***
*Ruonan Li, Ziyue Zhang, Xianjin Xia, et al. From Interference Mitigation to Toleration: Pathway to Practical Spatial Reuse in LPWANs, ACM MobiCom 2025*. 
## 📌 Hardware Requirements and Dependencies
### Hardware
- LoRa nodes: SX1276 + Arduino Uno (version: 1.8.X)

- Open source gateway: External + USRP N210 with more than 1 antenna (or other Software Defined Radios with equivalent functionality) + GNU Radio Platform

- Data trace processing: A workstation with Matlab (version: R2022a)
### Dependency
*This code is tested on Ubuntu 18.04*.
## HydraNet Spatial Reuse Protocol
Below is the protocol of HydraNet, illustrating how it achieve spatial reuse in LPWANs.

<div align="center">
    <img src="Spatial_reuse_protocol.jpg" alt="HydraNet Spatial Reuse Protocol" width="600">
</div>

## 📋 Folders
- **Lora_rf95_server:** We provide .ino (for LoRa COTS node to transmit and receive the probing packets and get the RSS).
  
- **LinkProbing:** We provide .grc (for GNU Radio to get the probing packets) and .m (for Matlab to estimate AoA and get the weight for each path).

- **Uplink** - provide the .grc (for GNU Radio to control the USRP to receive packets) and .m (for Matlab to demodulate) files for uplink opensource gateway process LoRa packets.

- **Downlink** - provide the .grc (for GNU Radio to control the USRP to transmit packets) and .m (for Matlab to generate power beam and combined packets .dat data) files for downlink opensource gateway process LoRa packets.
  
## 🛠 Deploying & Running HydraNet on a LoRa Testbed
### **To Run the Artifact with Provided Data Trace**  

- **Uplink:** We provide data traces in uplink that reproduce the key functions of HydraNet.

*Open the "Uplink" folder, and run exp_sec5_Hydranet_uplink_2x2.m, the detected packets and symbol error number are printed on the command window*.
- *Download the real-world uplink data (.dat) from the github release: https://github.com/LeonaRnLi/HydraNet/releases/tag/uplink*.
---
- **Downlink:** We provide data traces in downlink from Arduino monitor that reproduce how to process HydraNet downlink data and provide combined power beam packets generation functions. And the provided data are used for HydraNet spatial reuse performance evaluation (Fig. 14a one CDF line)
- *To evaluate real-world deployment PRR, open the **Downlink/Artifact/** folder and run exp_eval_plot_Spatial_Reuse_2_tx_CDF_PRR.m. This script computes the Packet Reception Ratio (PRR)*. ***In addition, we have provided an updated example plotting script, fig_eval_example_cdf_plot.m, which reproduces the CDF results as shown in Figure 14a of the paper. This script uses the real-world .mat files directly and produces the complete visualization, including formatting and legend layout.***.

***You may also modify this exp_eval_plot_Spatial_Reuse_2_tx_CDF_PRR.m file to compute the PRR of your own experiments, such as results from the COTS node monitor, for immediate evaluation***.

### **To Run the HydraNet on Your Testbed**
- **Uplink & Downlink** We provide .grc files and python framework for GNU Radio and .m for the power beamforming named example_classicMUSIC_RSS_2x2.m.
- **We give 2x2 Downlink deployment example** which is the main target scenario for HydraNet. And you need to download some downlink provided COTS nodes can recognized packets form the github release: https://github.com/LeonaRnLi/HydraNet/releases/tag/downlink
- ***Step:***  
*1. Open the "Lora_rf95_server" folder, and upload the server .ino on two LoRa SX1276 nodes for receive the downlink packets*.  <br>  
*2. Open the "Lora_rf95_server" folder, and upload the client .ino on LoRa nodes to simplify the link probing process, in this code, after upload the monitor can print the link attenuation*.<br>  
*3. After get the 4 link attenuation, use the example_classicMUSIC_RSS_2x2.m to get the GNU Radio Tx power*.<br>  
*4. Run the exp_sec5_beamforming.m with modified AoA angle parameter to get the combined packets*.<br>  
*5. Open the downlink/GNU_file folder, use beamforming_2x2.grc and verify the .dat file, the Tx power need to be set as Tx_power_2x2.json (you can enhace the Tx power with margin to make the system more robust)*.<br>  

## Video Demo
Coming soon...

## 📢 Contact & Citation
For questions, open an issue or contact:

- **Ruonan Li** - ruo-nan.li@connect.polyu.hk  
- **Ziyue Zhang** - ziyue.zhang@connect.polyu.hk  

If you use HydraNet in your research, please cite:

```bibtex
@inproceedings{hydranet2025,
  author    = {Ruonan Li and Ziyue Zhang and Xianjin Xia and Ningning Hou and 
               Wenchang Chai and Shiming Yu and Yuanqing Zheng and Tao Gu},
  title     = {From Interference Mitigation to Toleration: Pathway to Practical Spatial Reuse in LPWANs},
  booktitle = {Proceedings of the 31st Annual International Conference on Mobile Computing and Networking (ACM MobiCom ’25)},
  year      = {2025},
  publisher = {ACM},
  location  = {Hong Kong, China},
  doi       = {10.1145/3680207.3723483}
}
