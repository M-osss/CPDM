# Wiring

Signal is digital on both links: HDZero for video, ELRS CRSF for control. No 5.8 analog VTX.

```
                    XT60 12 AWG <= 60 mm
  6S LiHV 1050 150C --------------------> ESC pack pads + 1000 uF 50 V
                                              |
                                              +-- 8S-tolerant 5 V BEC --> FC
                                              +-- 10 V VTX pad ---------> HDZero 1W (use the 8–16 V pad, not 5 V)
                                              +-- motor phases ---------> XING2 2207 x4

  XR4 TX/RX/GND  -->  FC UART1  (serial_rx = CRSF, 420000)
  HDZero TX/RX    -->  FC UART2  (MSP DisplayPort + VTX)
  M10 GPS TX/RX   -->  FC UART6  (GPS, 115200)
  HDZero camera   -->  HDZero VTX CAM connector (stock coax, do not substitute analog cam)
```

## UART map (BLITZ Mini F7 / typical F7)

| UART | Protocol | Device |
|------|----------|--------|
| 1 | Serial RX, CRSF | XR4 Gemini Xrossband |
| 2 | MSP DisplayPort | HDZero 1W VTX |
| 6 | GPS | M10 |

Softserial is not used. HDZero needs a hardware UART.

## Antennas

- ELRS: two dual-band T-antennas. Left rear arm = 900 MHz path (XR4 ant 1). Right rear arm = 2.4 GHz path (XR4 ant 2). Polarization 90° relative. Carbon under the element kills the link; TPU standoff 15 mm off the arm.
- HDZero: single 5.8 GHz on the tail, 45° back, clear of the 2.4 element.
- Do not coil unused antenna pigtail. Trim nothing; use the stock dual-band Ts.

## Power notes

2550 KV on 6S will spike past 80 A pack on punch. XT30 is the wrong connector. XT60 with a short 12 AWG pigtail is the minimum. A 150 C 1050 mAh pack has a paper burst of 157 A; internal resistance, not the label, will sag first. If pack voltage under punch drops below ~21.5 V, the limiter is the battery, not the motors.

## Ground

Single-point pack negative. Do not ground the HDZero SMA shell to carbon.
