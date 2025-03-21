#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Classic Music Doa
# Generated: Wed Jul  3 17:50:13 2024
##################################################


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
from gnuradio import wxgui
from gnuradio.eng_option import eng_option
from gnuradio.fft import window
from gnuradio.filter import firdes
from gnuradio.filter import pfb
from gnuradio.wxgui import fftsink2
from gnuradio.wxgui import forms
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
import math
import time
import wx


class Classic_MUSIC_DOA(grc_wxgui.top_block_gui):

    def __init__(self):
        grc_wxgui.top_block_gui.__init__(self, title="Classic Music Doa")
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
        self.wxgui_fftsink2_0 = fftsink2.fft_sink_c(
        	self.GetWin(),
        	baseband_freq=0,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=samp_rate,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title='FFT Plot',
        	peak_hold=False,
        )
        self.Add(self.wxgui_fftsink2_0.win)
        self.uhd_usrp_source_ctrl_rx = uhd.usrp_source(
        	",".join(("addr0=192.168.10.4,addr1=192.168.10.5,addr2 =192.168.10.2,addr3=192.168.10.6", '')),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(4),
        	),
        )
        self.uhd_usrp_source_ctrl_rx.set_clock_source('external', 0)
        self.uhd_usrp_source_ctrl_rx.set_time_source('external', 0)
        self.uhd_usrp_source_ctrl_rx.set_clock_source('external', 1)
        self.uhd_usrp_source_ctrl_rx.set_time_source('external', 1)
        self.uhd_usrp_source_ctrl_rx.set_clock_source('external', 2)
        self.uhd_usrp_source_ctrl_rx.set_time_source('external', 2)
        self.uhd_usrp_source_ctrl_rx.set_clock_source('external', 3)
        self.uhd_usrp_source_ctrl_rx.set_time_source('external', 3)
        self.uhd_usrp_source_ctrl_rx.set_samp_rate(samp_rate)
        self.uhd_usrp_source_ctrl_rx.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_source_ctrl_rx.set_center_freq(frequency, 0)
        self.uhd_usrp_source_ctrl_rx.set_gain(rx_gain, 0)
        self.uhd_usrp_source_ctrl_rx.set_antenna('TX/RX', 0)
        self.uhd_usrp_source_ctrl_rx.set_center_freq(frequency, 1)
        self.uhd_usrp_source_ctrl_rx.set_gain(rx_gain, 1)
        self.uhd_usrp_source_ctrl_rx.set_antenna('TX/RX', 1)
        self.uhd_usrp_source_ctrl_rx.set_center_freq(frequency, 2)
        self.uhd_usrp_source_ctrl_rx.set_gain(rx_gain, 2)
        self.uhd_usrp_source_ctrl_rx.set_antenna('TX/RX', 2)
        self.uhd_usrp_source_ctrl_rx.set_center_freq(frequency, 3)
        self.uhd_usrp_source_ctrl_rx.set_gain(rx_gain, 3)
        self.uhd_usrp_source_ctrl_rx.set_antenna('TX/RX', 3)
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
        self.pfb_arb_resampler_zf_rx2 = pfb.arb_resampler_ccf(
        	  bw/samp_rate,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_zf_rx2.declare_sample_delay(0)

        self.pfb_arb_resampler_zf_rx1 = pfb.arb_resampler_ccf(
        	  bw/samp_rate,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_zf_rx1.declare_sample_delay(0)

        self.pfb_arb_resampler_ctrl_rx2 = pfb.arb_resampler_ccf(
        	  bw/samp_rate,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_ctrl_rx2.declare_sample_delay(0)

        self.pfb_arb_resampler_ctrl_rx1 = pfb.arb_resampler_ccf(
        	  bw/samp_rate,
                  taps=None,
        	  flt_size=32)
        self.pfb_arb_resampler_ctrl_rx1.declare_sample_delay(0)

        self.blocks_rotator_zf_rx2 = blocks.rotator_cc(0)
        self.blocks_rotator_zf_rx1 = blocks.rotator_cc(0)
        self.blocks_rotator_ce_rx2 = blocks.rotator_cc(0)
        self.blocks_rotator_ce_rx1 = blocks.rotator_cc(0)
        self.blocks_file_sink_zf_rx2 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment_Sec4_Design/RAWDATA/MUSIC_TEST/input/LoRa_MUSIC_RX4.dat', False)
        self.blocks_file_sink_zf_rx2.set_unbuffered(False)
        self.blocks_file_sink_zf_rx1 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment_Sec4_Design/RAWDATA/MUSIC_TEST/input/LoRa_MUSIC_RX3.dat', False)
        self.blocks_file_sink_zf_rx1.set_unbuffered(False)
        self.blocks_file_sink_ctrl_rx2 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment_Sec4_Design/RAWDATA/MUSIC_TEST/input/LoRa_MUSIC_RX2.dat', False)
        self.blocks_file_sink_ctrl_rx2.set_unbuffered(False)
        self.blocks_file_sink_ctrl_rx1 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/kevin/Ruonan_WorkSpace/20240405-RSS-CE-Experiment/Experiment_Sec4_Design/RAWDATA/MUSIC_TEST/input/LoRa_MUSIC_RX1.dat', False)
        self.blocks_file_sink_ctrl_rx1.set_unbuffered(False)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_rotator_ce_rx1, 0), (self.pfb_arb_resampler_ctrl_rx1, 0))
        self.connect((self.blocks_rotator_ce_rx2, 0), (self.pfb_arb_resampler_ctrl_rx2, 0))
        self.connect((self.blocks_rotator_zf_rx1, 0), (self.pfb_arb_resampler_zf_rx1, 0))
        self.connect((self.blocks_rotator_zf_rx2, 0), (self.pfb_arb_resampler_zf_rx2, 0))
        self.connect((self.pfb_arb_resampler_ctrl_rx1, 0), (self.blocks_file_sink_ctrl_rx1, 0))
        self.connect((self.pfb_arb_resampler_ctrl_rx1, 0), (self.wxgui_fftsink2_0, 0))
        self.connect((self.pfb_arb_resampler_ctrl_rx2, 0), (self.blocks_file_sink_ctrl_rx2, 0))
        self.connect((self.pfb_arb_resampler_zf_rx1, 0), (self.blocks_file_sink_zf_rx1, 0))
        self.connect((self.pfb_arb_resampler_zf_rx2, 0), (self.blocks_file_sink_zf_rx2, 0))
        self.connect((self.uhd_usrp_source_ctrl_rx, 0), (self.blocks_rotator_ce_rx1, 0))
        self.connect((self.uhd_usrp_source_ctrl_rx, 1), (self.blocks_rotator_ce_rx2, 0))
        self.connect((self.uhd_usrp_source_ctrl_rx, 2), (self.blocks_rotator_zf_rx1, 0))
        self.connect((self.uhd_usrp_source_ctrl_rx, 3), (self.blocks_rotator_zf_rx2, 0))

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
        self.wxgui_fftsink2_0.set_sample_rate(self.samp_rate)
        self.uhd_usrp_source_ctrl_rx.set_samp_rate(self.samp_rate)
        self.pfb_arb_resampler_zf_rx2.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_zf_rx1.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_ctrl_rx2.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_ctrl_rx1.set_rate(self.bw/self.samp_rate)

    def get_rx_gain(self):
        return self.rx_gain

    def set_rx_gain(self, rx_gain):
        self.rx_gain = rx_gain
        self._rx_gain_slider.set_value(self.rx_gain)
        self._rx_gain_text_box.set_value(self.rx_gain)
        self.uhd_usrp_source_ctrl_rx.set_gain(self.rx_gain, 0)

        self.uhd_usrp_source_ctrl_rx.set_gain(self.rx_gain, 1)

        self.uhd_usrp_source_ctrl_rx.set_gain(self.rx_gain, 2)

        self.uhd_usrp_source_ctrl_rx.set_gain(self.rx_gain, 3)


    def get_frequency(self):
        return self.frequency

    def set_frequency(self, frequency):
        self.frequency = frequency
        self._frequency_text_box.set_value(self.frequency)
        self.uhd_usrp_source_ctrl_rx.set_center_freq(self.frequency, 0)
        self.uhd_usrp_source_ctrl_rx.set_center_freq(self.frequency, 1)
        self.uhd_usrp_source_ctrl_rx.set_center_freq(self.frequency, 2)
        self.uhd_usrp_source_ctrl_rx.set_center_freq(self.frequency, 3)

    def get_bw(self):
        return self.bw

    def set_bw(self, bw):
        self.bw = bw
        self._bw_text_box.set_value(self.bw)
        self.pfb_arb_resampler_zf_rx2.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_zf_rx1.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_ctrl_rx2.set_rate(self.bw/self.samp_rate)
        self.pfb_arb_resampler_ctrl_rx1.set_rate(self.bw/self.samp_rate)


def main(top_block_cls=Classic_MUSIC_DOA, options=None):

    tb = top_block_cls()
    tb.Start(True)
    tb.Wait()


if __name__ == '__main__':
    main()
