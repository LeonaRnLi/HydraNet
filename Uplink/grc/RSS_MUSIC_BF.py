#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Rss Music Bf2
# Generated: Tue Jul 16 21:43:20 2024
##################################################

import ctypes
import sys
import subprocess
import logging
import time
import json
from lxml import etree
import matlab.engine
if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.filter import pfb
from gnuradio.wxgui import forms
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
import math
import time
import wx
MATLAB_CMD = '/usr/local/MATLAB/R2017a/bin/matlab'
SCRIPT1 = '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment-Sec5-Evaluation/exp_sec5_classicMUSIC_RSS_2x2.m'
FLAG_FILE = '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment-Sec5-Evaluation/Tx_power_2x2.json'
def run_matlab_script(matlab_script):
    logger.debug("Running Matlab script: {}".format(matlab_script))
    eng.run(matlab_script, nargout=0)
    logger.debug("Matlab script finished: {}".format(matlab_script))

def read_tx_power_from_json(json_file):
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)
            tx_power = data.get("TxPower", [0, 0])
            if len(tx_power) >= 2:
                tx_power1 = tx_power[0]
                tx_power2 = tx_power[1]
            else:
                tx_power1, tx_power2 = 0, 0
            return tx_power1, tx_power2
    except json.JSONDecodeError as e:
        logger.error("Error decoding JSON: {}".format(e))
    except FileNotFoundError:
        logger.error("JSON file not found: {}".format(json_file))
    except Exception as e:
        logger.error("Error reading JSON file: {}".format(e))
    return 0, 0


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
eng = matlab.engine.start_matlab()

class RSS_MUSIC_BF2(grc_wxgui.top_block_gui):

    def __init__(self):
        grc_wxgui.top_block_gui.__init__(self, title="Rss Music Bf")
        _icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
        self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

        ##################################################
        # Variables
        ##################################################
        self.tx_gain = tx_gain = 20
        self.samp_rate = samp_rate = 500e3
        self.rx_gain = rx_gain = 20
        self.frequency = frequency = 915e6
        self.bw = bw = 125e3

        ##################################################
        # Blocks
        ##################################################
        _rx_gain_sizer = wx.BoxSizer(wx.VERTICAL)
        self._rx_gain_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_rx_gain_sizer,
        	value=self.rx_gain,
        	callback=self.set_rx_gain,
        	label='rx_gain',
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._rx_gain_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_rx_gain_sizer,
        	value=self.rx_gain,
        	callback=self.set_rx_gain,
        	minimum=0,
        	maximum=89,
        	num_steps=90,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.Add(_rx_gain_sizer)
        self._frequency_text_box = forms.text_box(
        	parent=self.GetWin(),
        	value=self.frequency,
        	callback=self.set_frequency,
        	label='frequency',
        	converter=forms.float_converter(),
        )
        self.Add(self._frequency_text_box)
        self._bw_text_box = forms.text_box(
        	parent=self.GetWin(),
        	value=self.bw,
        	callback=self.set_bw,
        	label='bw',
        	converter=forms.float_converter(),
        )
        self.Add(self._bw_text_box)
        self.uhd_usrp_source_ctrl_rx2 = uhd.usrp_source(
        	",".join(("addr0=192.168.10.5,addr1=192.168.10.4", '')),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(2),
        	),
        )
        self.uhd_usrp_source_ctrl_rx2.set_clock_source('external', 0)
        self.uhd_usrp_source_ctrl_rx2.set_time_source('external', 0)
        self.uhd_usrp_source_ctrl_rx2.set_clock_source('external', 1)
        self.uhd_usrp_source_ctrl_rx2.set_time_source('external', 1)
        self.uhd_usrp_source_ctrl_rx2.set_samp_rate(samp_rate)
        self.uhd_usrp_source_ctrl_rx2.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_source_ctrl_rx2.set_center_freq(frequency, 0)
        self.uhd_usrp_source_ctrl_rx2.set_gain(rx_gain, 0)
        self.uhd_usrp_source_ctrl_rx2.set_antenna('TX/RX', 0)
        self.uhd_usrp_source_ctrl_rx2.set_center_freq(frequency, 1)
        self.uhd_usrp_source_ctrl_rx2.set_gain(rx_gain, 1)
        self.uhd_usrp_source_ctrl_rx2.set_antenna('TX/RX', 1)
        self.uhd_usrp_source_ctrl_rx1 = uhd.usrp_source(
        	",".join(("addr0=192.168.10.6,addr1=192.168.10.2", '')),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(2),
        	),
        )
        self.uhd_usrp_source_ctrl_rx1.set_clock_source('external', 0)
        self.uhd_usrp_source_ctrl_rx1.set_time_source('external', 0)
        self.uhd_usrp_source_ctrl_rx1.set_clock_source('external', 1)
        self.uhd_usrp_source_ctrl_rx1.set_time_source('external', 1)
        self.uhd_usrp_source_ctrl_rx1.set_samp_rate(samp_rate)
        self.uhd_usrp_source_ctrl_rx1.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_source_ctrl_rx1.set_center_freq(frequency, 0)
        self.uhd_usrp_source_ctrl_rx1.set_gain(rx_gain, 0)
        self.uhd_usrp_source_ctrl_rx1.set_antenna('TX/RX', 0)
        self.uhd_usrp_source_ctrl_rx1.set_center_freq(frequency, 1)
        self.uhd_usrp_source_ctrl_rx1.set_gain(rx_gain, 1)
        self.uhd_usrp_source_ctrl_rx1.set_antenna('TX/RX', 1)
        self.uhd_usrp_sink_BFTx2 = uhd.usrp_sink(
        	",".join(("addr0=192.168.10.5,addr1=192.168.10.4", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(2),
        	),
        )
        self.uhd_usrp_sink_BFTx2.set_clock_source('external', 0)
        self.uhd_usrp_sink_BFTx2.set_time_source('external', 0)
        self.uhd_usrp_sink_BFTx2.set_clock_source('external', 1)
        self.uhd_usrp_sink_BFTx2.set_time_source('external', 1)
        self.uhd_usrp_sink_BFTx2.set_samp_rate(samp_rate)
        self.uhd_usrp_sink_BFTx2.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_sink_BFTx2.set_center_freq(frequency, 0)
        self.uhd_usrp_sink_BFTx2.set_gain(tx_gain, 0)
        self.uhd_usrp_sink_BFTx2.set_center_freq(frequency, 1)
        self.uhd_usrp_sink_BFTx2.set_gain(tx_gain, 1)
        self.uhd_usrp_sink_BFTx1 = uhd.usrp_sink(
        	",".join(("addr0=192.168.10.6,addr1=192.168.10.2", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(2),
        	),
        )
        self.uhd_usrp_sink_BFTx1.set_clock_source('external', 0)
        self.uhd_usrp_sink_BFTx1.set_time_source('external', 0)
        self.uhd_usrp_sink_BFTx1.set_clock_source('external', 1)
        self.uhd_usrp_sink_BFTx1.set_time_source('external', 1)
        self.uhd_usrp_sink_BFTx1.set_samp_rate(samp_rate)
        self.uhd_usrp_sink_BFTx1.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_sink_BFTx1.set_center_freq(frequency, 0)
        self.uhd_usrp_sink_BFTx1.set_gain(tx_gain, 0)
        self.uhd_usrp_sink_BFTx1.set_center_freq(frequency, 1)
        self.uhd_usrp_sink_BFTx1.set_gain(tx_gain, 1)
        _tx_gain_sizer = wx.BoxSizer(wx.VERTICAL)
        self._tx_gain_text_box = forms.text_box(
        	parent=self.GetWin(),
        	sizer=_tx_gain_sizer,
        	value=self.tx_gain,
        	callback=self.set_tx_gain,
        	label='tx_gain',
        	converter=forms.float_converter(),
        	proportion=0,
        )
        self._tx_gain_slider = forms.slider(
        	parent=self.GetWin(),
        	sizer=_tx_gain_sizer,
        	value=self.tx_gain,
        	callback=self.set_tx_gain,
        	minimum=0,
        	maximum=89,
        	num_steps=90,
        	style=wx.SL_HORIZONTAL,
        	cast=float,
        	proportion=1,
        )
        self.Add(_tx_gain_sizer)
        self.pfb_arb_resampler_tx22 = pfb.arb_resampler_ccf(
        	  samp_rate/bw,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_tx22.declare_sample_delay(0)

        self.pfb_arb_resampler_tx21 = pfb.arb_resampler_ccf(
        	  samp_rate/bw,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_tx21.declare_sample_delay(0)

        self.pfb_arb_resampler_tx12 = pfb.arb_resampler_ccf(
        	  samp_rate/bw,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_tx12.declare_sample_delay(0)

        self.pfb_arb_resampler_tx11 = pfb.arb_resampler_ccf(
        	  samp_rate/bw,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_tx11.declare_sample_delay(0)

        self.pfb_arb_resampler_rx22 = pfb.arb_resampler_ccf(
        	  bw/samp_rate,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_rx22.declare_sample_delay(0)

        self.pfb_arb_resampler_rx21 = pfb.arb_resampler_ccf(
        	  bw/samp_rate,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_rx21.declare_sample_delay(0)

        self.pfb_arb_resampler_rx12 = pfb.arb_resampler_ccf(
        	  bw/samp_rate,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_rx12.declare_sample_delay(0)

        self.pfb_arb_resampler_rx11 = pfb.arb_resampler_ccf(
        	  bw/samp_rate,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_rx11.declare_sample_delay(0)


        self.blocks_rotator_rx22 = blocks.rotator_cc(0)
        self.blocks_rotator_rx21 = blocks.rotator_cc(0)
        self.blocks_rotator_rx11 = blocks.rotator_cc(0)
        self.blocks_rotator_rx12 = blocks.rotator_cc(0)
       
        self.blocks_file_sink_rx22 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment-Sec5-Evaluation/RAWDATA/MUSIC_TEST/input/LoRa_MUSIC_RX4.dat', False)
        self.blocks_file_sink_rx22.set_unbuffered(False)
        self.blocks_file_sink_rx21 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment-Sec5-Evaluation/RAWDATA/MUSIC_TEST/input/LoRa_MUSIC_RX3.dat', False)
        self.blocks_file_sink_rx21.set_unbuffered(False)
        self.blocks_file_sink_rx12 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment-Sec5-Evaluation/RAWDATA/MUSIC_TEST/input/LoRa_MUSIC_RX2.dat', False)
        self.blocks_file_sink_rx12.set_unbuffered(False)
        self.blocks_file_sink_rx11 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment-Sec5-Evaluation/RAWDATA/MUSIC_TEST/input/LoRa_MUSIC_RX1.dat', False)
        self.blocks_file_sink_rx11.set_unbuffered(False)

       
                ##################################################
        # Connections
        ##################################################
        self.connect((self.uhd_usrp_source_ctrl_rx1, 0), (self.pfb_arb_resampler_rx11, 0))
        self.connect((self.uhd_usrp_source_ctrl_rx1, 1), (self.pfb_arb_resampler_rx12, 0))
        self.connect((self.uhd_usrp_source_ctrl_rx2, 0), (self.pfb_arb_resampler_rx21, 0))
        self.connect((self.uhd_usrp_source_ctrl_rx2, 1), (self.pfb_arb_resampler_rx22, 0))

        ##################################################
        # Iterations
        ##################################################
        self.blocks_rotator_tx11 = blocks.rotator_cc(0)
        self.blocks_rotator_tx12 = blocks.rotator_cc(0)
        self.blocks_rotator_tx21 = blocks.rotator_cc(0)
        self.blocks_rotator_tx22 = blocks.rotator_cc(0)

        self.blocks_file_source_BF22 = blocks.file_source(gr.sizeof_gr_complex * 1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment-Sec5-Evaluation/RAWDATA/LoRa_CtrlPKT/Tx_packet/Beamforming_A2LTx2.dat', True)
        self.blocks_file_source_BF21 = blocks.file_source(gr.sizeof_gr_complex * 1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment-Sec5-Evaluation/RAWDATA/LoRa_CtrlPKT/Tx_packet/Beamforming_A2LTx1.dat', True)
        self.blocks_file_source_BF12 = blocks.file_source(gr.sizeof_gr_complex * 1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment-Sec5-Evaluation/RAWDATA/LoRa_CtrlPKT/Tx_packet/Beamforming_1to9Tx2.dat', True)
        self.blocks_file_source_BF11 = blocks.file_source(gr.sizeof_gr_complex * 1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment-Sec5-Evaluation/RAWDATA/LoRa_CtrlPKT/Tx_packet/Beamforming_1to9Tx1.dat', True)

        self.connect((self.blocks_file_source_BF11, 0), (self.blocks_rotator_tx11, 0))
        self.connect((self.blocks_file_source_BF12, 0), (self.blocks_rotator_tx12, 0))
        self.connect((self.blocks_file_source_BF21, 0), (self.blocks_rotator_tx21, 0))
        self.connect((self.blocks_file_source_BF22, 0), (self.blocks_rotator_tx22, 0))
        self.connect((self.blocks_rotator_tx11, 0), (self.pfb_arb_resampler_tx11, 0))
        self.connect((self.blocks_rotator_tx12, 0), (self.pfb_arb_resampler_tx12, 0))
        self.connect((self.blocks_rotator_tx21, 0), (self.pfb_arb_resampler_tx21, 0))
        self.connect((self.blocks_rotator_tx22, 0), (self.pfb_arb_resampler_tx22, 0))
        self.connect((self.pfb_arb_resampler_tx11, 0), (self.uhd_usrp_sink_BFTx1, 0))
        self.connect((self.pfb_arb_resampler_tx12, 0), (self.uhd_usrp_sink_BFTx1, 1))
        self.connect((self.pfb_arb_resampler_tx21, 0), (self.uhd_usrp_sink_BFTx2, 0))
        self.connect((self.pfb_arb_resampler_tx22, 0), (self.uhd_usrp_sink_BFTx2, 1))

        def run_iteration():
            self.start()
            time.sleep(2)  # 运行时间
            self.stop()
            self.wait()

            logger.debug("waiting for the MATLAB Tx_power Value...")
            run_matlab_script(SCRIPT1)

            # 从JSON文件读取发射功率值
            tx_power1, tx_power2 = read_tx_power_from_json(FLAG_FILE)
            logger.info("Tx Power1: {}, Tx Power2: {}".format(tx_power1, tx_power2))

            # 更新发射功率
            self.uhd_usrp_sink_BFTx1.set_gain(tx_power1, 0)
            self.uhd_usrp_sink_BFTx1.set_gain(tx_power1, 1)
            self.uhd_usrp_sink_BFTx2.set_gain(tx_power2, 0)
            self.uhd_usrp_sink_BFTx2.set_gain(tx_power2, 1)

            self.start()
            time.sleep(0.3)  # 运行时间
            self.stop()
            self.wait()

        for i in range(3):  # 循环三次
            logger.info("Iteration: {}".format(i + 1))
            run_iteration()

        # 停止脚本和 USRP
        self.stop()
        self.wait()

 

    def get_tx_gain(self):
        return self.tx_gain

    def set_tx_gain(self, tx_gain):
        self.tx_gain = tx_gain
        self._tx_gain_slider.set_value(self.tx_gain)
        self._tx_gain_text_box.set_value(self.tx_gain)

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_source_ctrl_rx2.set_samp_rate(self.samp_rate)
        self.uhd_usrp_source_ctrl_rx1.set_samp_rate(self.samp_rate)
        self.uhd_usrp_sink_BFTx2.set_samp_rate(self.samp_rate)
        self.uhd_usrp_sink_BFTx1.set_samp_rate(self.samp_rate)
        self.pfb_arb_resampler_tx22.set_rate(self.samp_rate/self.bw)
        self.pfb_arb_resampler_tx21.set_rate(self.samp_rate/self.bw)
        self.pfb_arb_resampler_tx12.set_rate(self.samp_rate/self.bw)
        self.pfb_arb_resampler_tx11.set_rate(self.samp_rate/self.bw)
        self.pfb_arb_resampler_rx22.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_rx21.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_rx12.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_rx11.set_rate(self.bw/self.samp_rate)

    def get_rx_gain(self):
        return self.rx_gain

    def set_rx_gain(self, rx_gain):
        self.rx_gain = rx_gain
        self._rx_gain_slider.set_value(self.rx_gain)
        self._rx_gain_text_box.set_value(self.rx_gain)
        self.uhd_usrp_source_ctrl_rx2.set_gain(self.rx_gain, 0)

        self.uhd_usrp_source_ctrl_rx2.set_gain(self.rx_gain, 1)

        self.uhd_usrp_source_ctrl_rx2.set_gain(self.rx_gain, 2)

        self.uhd_usrp_source_ctrl_rx2.set_gain(self.rx_gain, 3)

        self.uhd_usrp_source_ctrl_rx1.set_gain(self.rx_gain, 0)

        self.uhd_usrp_source_ctrl_rx1.set_gain(self.rx_gain, 1)

        self.uhd_usrp_source_ctrl_rx1.set_gain(self.rx_gain, 2)

        self.uhd_usrp_source_ctrl_rx1.set_gain(self.rx_gain, 3)


    def get_frequency(self):
        return self.frequency

    def set_frequency(self, frequency):
        self.frequency = frequency
        self._frequency_text_box.set_value(self.frequency)
        self.uhd_usrp_source_ctrl_rx2.set_center_freq(self.frequency, 0)
        self.uhd_usrp_source_ctrl_rx2.set_center_freq(self.frequency, 1)
        self.uhd_usrp_source_ctrl_rx2.set_center_freq(self.frequency, 2)
        self.uhd_usrp_source_ctrl_rx2.set_center_freq(self.frequency, 3)
        self.uhd_usrp_source_ctrl_rx1.set_center_freq(self.frequency, 0)
        self.uhd_usrp_source_ctrl_rx1.set_center_freq(self.frequency, 1)
        self.uhd_usrp_source_ctrl_rx1.set_center_freq(self.frequency, 2)
        self.uhd_usrp_source_ctrl_rx1.set_center_freq(self.frequency, 3)

    def get_bw(self):
        return self.bw

    def set_bw(self, bw):
        self.bw = bw
        self._bw_text_box.set_value(self.bw)
        self.pfb_arb_resampler_tx22.set_rate(self.samp_rate/self.bw)
        self.pfb_arb_resampler_tx21.set_rate(self.samp_rate/self.bw)
        self.pfb_arb_resampler_tx12.set_rate(self.samp_rate/self.bw)
        self.pfb_arb_resampler_tx11.set_rate(self.samp_rate/self.bw)
        self.pfb_arb_resampler_rx22.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_rx21.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_rx12.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_rx11.set_rate(self.bw/self.samp_rate)


def main(top_block_cls=RSS_MUSIC_BF2, options=None):

    tb = top_block_cls()
    tb.Start(True)
    tb.Wait()


if __name__ == '__main__':
    main()
