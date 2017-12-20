# Read common PicoHarp data files and optionally write them to disk.
# Based on the Original matlab and C code by
# Peter Kapusta, Michael Wahl, PicoQuant GmbH 2006, updated May 2007
# Converted to python by Bernd Sontheimer, Humboldt University Berlin, August 2017
# Use at your own risk. No warranties.
# Data is kept in memory despite of the huge amount of memory
# this can take in case of large files. Therefore the file can be split using
# the records_per_split parameter.
# Note that marker events have a lower time resolution and may therefore appear
# in the file slightly out of order with respect to regular (photon) event records.
# This is by design. Markers are designed only for relatively coarse
# synchronization requirements such as image scanning.

import numpy as _np
import time


def read_picoquant_legacy_header(fid):
    """
    read a picoharp legacy data file header.
    :param fid: input file handle
    :return: file header
    """
    ascii_header = dict()
    binary_header = dict()
    board_header = dict()
    t_mode_header = dict()

    # read ascii header
    ascii_header['Identifier'] = fid.read(16).decode('utf-8').replace('\x00', '')
    ascii_header['Format Version'] = fid.read(6).decode('utf-8').replace('\x00', '')
    if ascii_header['Format Version'] != '2.0':
        print('Warning: This program is for version 2.0 files only.')
    ascii_header['Creator Name'] = fid.read(18).decode('utf-8').replace('\x00', '')
    ascii_header['Creator Version'] = fid.read(12).decode('utf-8').replace('\x00', '')
    ascii_header['File Time'] = fid.read(18).decode('utf-8').replace('\x00', '')
    ascii_header['CRLF'] = fid.read(2).decode('utf-8').replace('\x00', '')
    ascii_header['Comment'] = fid.read(256).decode('utf-8').replace('\x00', '')

    # read binary header
    binary_header['Curves'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Bits per Record'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Rounting Channels'], =  _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Number of Boards'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Active Curve'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Measurement Mode'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Sub-Mode'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Range No.'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Offset / ns'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Acquisition Time / ms'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Stop at'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Stop on Overflow'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Restart'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Lin/Log'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Time Axis from'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Time Axis to'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Count Axis from'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Count Axis to'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Curve map to'] = _np.fromfile(fid, count=8, dtype=_np.int32)
    binary_header['Display Curve show'] = _np.fromfile(fid, count=8, dtype=_np.int32)
    binary_header['Param0'] = {'Start': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'Stop': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'End': _np.fromfile(fid, count=1, dtype=_np.float32)[0]
                               }
    binary_header['Param1'] = {'Start': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'Stop': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'End': _np.fromfile(fid, count=1, dtype=_np.float32)[0]
                               }
    binary_header['Param2'] = {'Start': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'Stop': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'End': _np.fromfile(fid, count=1, dtype=_np.float32)[0]
                               }
    binary_header['Repeat Mode'], =  _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Repeats per Curve'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Repeat Time'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Repeat Wait Time'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Script Name'] = fid.read(20).decode('utf-8').replace('\x00', '')

    # read board specific header
    board_header['Hardware Identifier'] = fid.read(16).decode('utf-8').replace('\x00', '')
    board_header['Hardware Version'] = fid.read(8).decode('utf-8').replace('\x00', '')
    board_header['Hardware Serial'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['Sync Divider'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['CFD ZeroCross (Ch0) / mV'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['CFD Dischr. (Ch0) / mV'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['CFD ZeroCross (Ch1) / mV'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['CFD Dischr. (Ch1) / mV'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['Resolution / ns'], = _np.fromfile(fid, count=1, dtype=_np.float32)
    board_header['Router Model Code'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['Router Enabled'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['Rounter Ch1'] = {'Input Type': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Edge': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Present': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD ZeroCross / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   }
    board_header['Rounter Ch2'] = {'Input Type': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Edge': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Present': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD ZeroCross / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   }
    board_header['Rounter Ch3'] = {'Input Type': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Edge': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Present': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD ZeroCross / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   }
    board_header['Rounter Ch4'] = {'Input Type': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Edge': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Present': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD ZeroCross / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   }

    # Router settings are meaningful only for an existing router:
    if board_header['Router Model Code'] == 0:
        del board_header['Rounter Ch1'],
        del board_header['Rounter Ch2'],
        del board_header['Rounter Ch3'],
        del board_header['Rounter Ch4']

    # Read T mode specific header
    t_mode_header['External Devices'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Reserved1'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Reserved2'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Count Rate (Ch0) / Hz'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Count Rate (Ch1) / Hz'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Stop after / ms'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Stop Reason'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Number of Records'], = _np.fromfile(fid, count=1, dtype=_np.uint32)

    # Read special imaging header
    imaging_header_size, = _np.fromfile(fid, count=1, dtype=_np.uint32)
    img_header = _np.fromfile(fid, count=imaging_header_size, dtype=_np.int32)

    # combine headers
    header = {
        'ASCII': ascii_header,
        'Binary': binary_header,
        'Board': board_header,
        'T_mode': t_mode_header,
        'Image': img_header
    }
    return header


def read_ptu_header(fid):
    """
    read a picoharp .ptu data file header.
    :param fid: input file handle
    :return: file header
    """
    # Tag Types
    tag_types = {
        0xFFFF0008: 'tyEmpty8',
        0x00000008: 'tyBool8',
        0x10000008: 'tyInt8',
        0x11000008: 'tyBitSet64',
        0x12000008: 'tyColor8',
        0x20000008: 'tyFloat8',
        0x21000008: 'tyTDateTime',
        0x2001FFFF: 'tyFloat8Array',
        0x4001FFFF: 'tyAnsiString',
        0x4002FFFF: 'tyWideString',
        0xFFFFFFFF: 'tyBinaryBlob'
    }

    # Record Types
    rec_types = {
        0x00010303: 'rtPicoHarpT3',
        0x00010203: 'rtPicoHarpT2',
        0x00010304: 'rtHydraHarpT3',
        0x00010204: 'rtHydraHarpT2',
        0x01010304: 'rtHydraHarp2T3',
        0x01010204: 'rtHydraHarp2T2',
        0x00010305: 'rtTimeHarp260NT3',
        0x00010205: 'rtTimeHarp260NT2',
        0x00010306: 'rtTimeHarp260PT3',
        0x00010206: 'rtTimeHarp260PT2'
    }

    header = dict()
    magic = fid.read(8).decode('utf-8').replace('\x00', '')
    if magic != 'PQTTTR':
        raise IOError('Not a valid PTU file. '
                      'Magic: '.format(magic))
    header['Magic'] = magic
    header['Version'] = fid.read(8).decode('utf-8').replace('\x00', '')

    while True:
        # read tag header
        tag_ident = fid.read(32).decode('utf-8').replace('\x00', '')
        tag_idx = _np.fromfile(fid, count=1, dtype=_np.int32)[0]
        tag_type = tag_types[_np.fromfile(fid, count=1, dtype=_np.uint32)[0]]
        tag_data = None

        # most tag_value types are int64 but some are double
        if (tag_type == 'tyFloat8') | (type == 'tyFloat8Array'):
            dtype = _np.double
        else:
            dtype = _np.int64

        tag_value = _np.fromfile(fid, count=1, dtype=dtype)[0]

        # Some tag types need additional conversion
        if tag_type == 'tyFloat8':
            tag_value = _np.int64(tag_value).view(_np.float64)
        elif tag_type == 'tyBool8':
            tag_value = bool(tag_value)
        elif tag_type == 'tyTDateTime':
            tag_value = convert_ptu_time(_np.uint64(tag_value).view(_np.float64))

        # Some tag types have additional tag_data
        if tag_type == 'tyAnsiString':
            tag_data = fid.read(tag_value).decode('utf-8').replace('\x00', '')
        elif tag_type == 'tyFloat8Array':
            tag_data = fid.read(tag_value).decode('utf-8').replace('\x00', '')
        elif tag_type == 'tyWideString':
            tag_data = fid.read(tag_value).decode('utf-8').replace('\x00', '')
        elif tag_type == 'tyBinaryBlob':
            tag_data = fid.read(tag_value)
        if tag_ident == "Header_End":
            break
        header[tag_ident] = dict(tag_idx=tag_idx, tag_value=tag_value, tag_data=tag_data)

    # convert to readable rec type
    rec_type = rec_types[header['TTResultFormat_TTTRRecType']['tag_value']]
    header['TTResultFormat_TTTRRecType']['tag_value'] = rec_type
    return header


def read_pt2_records(fid, records):
    """
    read a picoharp T2 Mode records file header.
    :param fid: input file handle
    :param records: number of records to read
    :return: data: record channel and timestamp numpy arrays
    """
    # Read the T2 mode event records
    resolution = 4e-12  # 4ps is the standard max resolution of the picoharp
    wraparound = _np.uint64(210698240)
    t2_records = _np.fromfile(fid, count=records, dtype=_np.uint32)  # each record is composed of 32 bits
    t2_time = _np.bitwise_and(268435455, t2_records)  # last 28 bits are the time stamp
    t2_channel = _np.right_shift(t2_records, 28)  # first 4 bits are the channel
    t2_channel = _np.uint8(t2_channel)  # first 4 bits are the channel
    del t2_records
    # if the channel is 15, then the record is special
    marker_mask = t2_channel == 15
    # For special records, the lowest 4 bits of the time stamp are actually marker bits
    # if this is 0, the marker is an overflow marker
    ofl_marker = marker_mask & (_np.bitwise_and(15, t2_time) == 0)
    t2_time = _np.uint64(t2_time)
    t2_time += _np.cumsum(ofl_marker, dtype=_np.uint64) * wraparound
    t2_time = t2_time * resolution
    # remove overflows from data:
    t2_channel = t2_channel[~ofl_marker]
    t2_time = t2_time[~ofl_marker]
    return [t2_channel, t2_time]


def read_pt3_records(fid, records):
    """
    read a picoharp T3 Mode records file header.
    :param fid: input file handle
    :param records: number of records to read
    :return: data: record channel and timestamp numpy arrays
    """
    resolution = 4e-12  # 4 ps is the standard max resolution of the picoharp
    # wraparound = _np.uint64(65536)
    t3_records = _np.fromfile(fid, count=records, dtype=_np.uint32)  # each record is composed of 32 bits
    # for now the sync time stamp is ignored
    # t3_sync = _np.bitwise_and(65535, t3_records)  # last 16 bits are the sync counter time stamp
    t3_time = _np.bitwise_and(_np.right_shift(t3_records, 16), 4095)
    t3_channel = _np.right_shift(t3_records, 28)  # first 4 bits are the channel
    t3_channel = _np.uint8(t3_channel)
    del t3_records  # free up unused memory
    # if the channel is 15, then the record is special
    marker_mask = t3_channel == 15
    # For special records, the lowest 4 bits of the time stamp are actually marker bits
    # if this is 0, the marker is an overflow marker
    ofl_marker = marker_mask & (_np.bitwise_and(15, t3_time) == 0)
    # remove overflows from data:
    t3_channel = t3_channel[~ofl_marker]
    t3_time = t3_time[~ofl_marker]
    t3_time = t3_time * resolution
    # ofl_correction = _np.cumsum(ofl_marker, dtype=_np.uint64) * wraparound
    # t3_sync = _np.uint64(t3_sync)
    # t3_sync += ofl_correction
    return [t3_channel, t3_time]


def read_picoquant(s, records_per_split=_np.infty, save_as_npz=False):
    """
    read a picoharp data files. As of now supported are .pt2 .pt3 and some .ptu files.
    :param s: input file string.
    :param records_per_split: big files may be split if they do not fit into memory.
    if set, save_as_npz defaults to True.
    :param save_as_npz: default is False. if set True the data is written to .npz files in the same directory as a the
    input file.
    :return header: dict of the file header information
    :return data: data converted to list of numpy arrays. Entries are dependent on the input file
    """
    with open(s, 'rb') as fid:
        file_ending = s[-4:]
        # read headers and number of records
        if (file_ending == '.pt2') | (file_ending == '.pt3'):
            header = read_picoquant_legacy_header(fid)
            number_of_records = header['T_mode']['Number of Records']

        elif file_ending == '.ptu':
            header = read_ptu_header(fid)
            number_of_records = header['TTResult_NumberOfRecords']['tag_value']
        else:
            print('not a valid file. Only .pt2, .pt3 and .ptu files are supported!')
            return False

        if number_of_records <= records_per_split:
            records_to_read = [number_of_records]
        else:
            splits, remainder = divmod(number_of_records, records_per_split)
            records_to_read = [records_per_split]*splits + [remainder]
            save_as_npz = True
            print('Save output and split file into {} parts'.format(splits))

        data = None
        for i, records in enumerate(records_to_read):
            if file_ending == '.pt2':
                data = read_pt2_records(fid, records)
            if file_ending == '.pt3':
                data = read_pt3_records(fid, records)
            if file_ending == '.ptu':
                rec_type = header['TTResultFormat_TTTRRecType']['tag_value']
                if rec_type == 'rtPicoHarpT3':
                    data = read_pt3_records(fid, records)
                elif rec_type == 'rtPicoHarpT2':
                    data = read_pt2_records(fid, records)
                else:
                    print('Illegal or not implemented RecordType')
                    return False
            if save_as_npz:
                _np.savez(s[:-4] + str(i), header=header, data=data)
        return header, data


def convert_ptu_time(tdatetime):
    """Convert the time used in PTU files to nice time string."""
    epoch_diff = 25569  # days between 30/12/1899 and 01/01/1970
    secs_in_day = 86400  # number of seconds in a day
    t = time.gmtime((tdatetime - epoch_diff) * secs_in_day)
    t = time.strftime("%Y/%m/%d %H:%M:%S", t)
    return t
