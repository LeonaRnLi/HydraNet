// rf95_client.pde
// -*- mode: C++ -*-
// Example sketch showing how to create a simple messageing client
// with the RH_RF95 class. RH_RF95 class does not provide for addressing or
// reliability, so you should only use RH_RF95 if you do not need the higher
// level messaging abilities.
// It is designed to work with the other example rf95_server
// Tested with Anarduino MiniWirelessLoRa, Rocket Scream Mini Ultra Pro with
// the RFM95W, Adafruit Feather M0 with RFM95

#include <SPI.h>
#include <RH_RF95.h>

#define ISADAFRUIT 0 //0 for shield, 1 for adafruit

#if ISADAFRUIT == 1 
  RH_RF95 rf95(8, 3); // Adafruit Feather M0 with RFM95 
#else
  RH_RF95 rf95;  
#endif

// Singleton instance of the radio driver
//RH_RF95 rf95;
//RH_RF95 rf95(5, 2); // Rocket Scream Mini Ultra Pro with the RFM95W
//RH_RF95 rf95(8, 3); // Adafruit Feather M0 with RFM95 

// Need this on Arduino Zero with SerialUSB port (eg RocketScream Mini Ultra Pro)
//#define Serial SerialUSB

static uint8_t count = 0;
int led = LED_BUILTIN;
void setup() 
{
  // Rocket Scream Mini Ultra Pro with the RFM95W only:
  // Ensure serial flash is not interfering with radio communication on SPI bus
  pinMode(led, OUTPUT);
   
  Serial.begin(9600);
  
  //while (!Serial) ; // Wait for serial port to be available
  
  if (!rf95.init())
    Serial.println("init failed");

  // Defaults after init are 434.0MHz, 13dBm, Bw = 125 kHz, Cr = 4/5, Sf = 128chips/symbol, CRC on
  
  // The default transmitter power is 13dBm, using PA_BOOST.
  // If you are using RFM95/96/97/98 modules which uses the PA_BOOST transmitter pin, then 
  // you can set transmitter powers from 5 to 23 dBm:
  //  driver.setTxPower(23, false);
  // If you are using Modtronix inAir4 or inAir9,or any other module which uses the
  // transmitter RFO pins and not the PA_BOOST pins
  // then you can configure the power transmitter power for -1 to 14 dBm and with useRFO true. 
  // Failure to do that will result in extremely low transmit powers.
  //  driver.setTxPower(14, true);

  // Config LoRa Radio
  // Setup ISM frequency
  rf95.setFrequency(915);
    ////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Setup Spreading Factor (6 ~ 12)
  rf95.setSpreadingFactor(8);


  // Setup BandWidth, option: 7800,10400,15600,20800,31250,41700,62500,125000,250000,500000
  //Lower BandWidth for longer distance.
  rf95.setSignalBandwidth(250E3); // Should be same in USRP and LORA SERVER, 125e3 also works if set same
   ////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
  // Setup Coding Rate:5(4/5),6(4/6),7(4/7),8(4/8) 
  rf95.setCodingRate4(5);
  
  // Setup Power,dBm
  rf95.setTxPower(23,false);
  rf95.setPreambleLength(8);
  rf95.setImplicitHeaderMode(true); // Matt's gr-lora only support Implicite Mode: The whitening sequence for Implicit and Explicte modes are different
  rf95.setLowDataRateMobileNodeMode(true); // Matt's gr-lora only support LDR
  rf95.setEnableCRC(false); // USRP decodes CRC, Should be set same with LORA SERVER (CRC on/off)
  
  rf95.printRegisters();
 
  Serial.println("Waiting for radio to setup");
  delay(1000);
  Serial.println("Client setup completed");

    rf95.printRegisters();

}

void loop()
{
  #if ISADAFRUIT == 1 
    uint8_t data[] = "Ada-Client";
  #else
    uint8_t data[] = "123456789#";//ABCDEFGHJKL
    //uint8_t data[] = {0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x00};
  #endif     
  
  //uint8_t data[] = {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
  //0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00};

  //uint8_t data[] = {0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,
  //0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff};
  
  char buffer[RH_RF95_MAX_MESSAGE_LEN];
  //sprintf(buffer, "%d: %s\0", count++, data);
  sprintf(buffer, "%s\0", data);
  
  // Send a message to rf95_server
  #if ISADAFRUIT == 1 
    Serial.print("AdaClient sending - ");
  #else
    Serial.print("Shield sending - ");
  #endif 
  
  Serial.println(buffer);
  Serial.println(strlen(buffer));
  //rf95.send((uint8_t*)buffer, sizeof(buffer));
  rf95.send((uint8_t*)buffer, strlen(buffer)+1);
  rf95.waitPacketSent();
  //blink_send(led);

  ///*
  // Now wait for a reply
  uint8_t buf[RH_RF95_MAX_MESSAGE_LEN];
  uint8_t len = sizeof(buf);

  if (rf95.waitAvailableTimeout(4000))
  { 
    // Should be a reply message for us now   
    if (rf95.recv(buf, &len))
   {
      Serial.print("Client got packet from SERVER: ");
      Serial.println((char*)buf);
      blink_receive(led);
      Serial.print("RSSI: ");
      Serial.println(rf95.lastRssi(), DEC);    
    }
    else
    {
      Serial.println("recv failed");
    }
  }
  else
  {
    // Serial.println("No reply from SERVER");
  }
  delay(3000/5);
  //*/
}

void blink_send(int led){
  digitalWrite(led, HIGH);
  delay(500);
  digitalWrite(led, LOW);
  delay(500);  
}

void blink_receive(int led){
  digitalWrite(led, HIGH);
  delay(250);
  digitalWrite(led, LOW);
  delay(250);
  digitalWrite(led, HIGH);
  delay(250);
  digitalWrite(led, LOW);
  delay(250);
}

/*
static uint8_t WhiteningKeyMSB= 0x01; // Global variable so the value is kept after starting the
static uint8_t WhiteningKeyLSB= 0xFF; // de-whitening process
void Sx1272RadioComputeWhitening( uint8_t *buffer, uint16_t bufferSize )
{
 uint8_t i;
 uint16_t j;
 uint8_t WhiteningKeyMSBPrevious;

 j = 0;
 buffer[j] ^= WhiteningKeyLSB;
 for( i = 0; i < 9; i++ )
 {
 WhiteningKeyMSBPrevious = WhiteningKeyMSB;
 WhiteningKeyMSB = ( WhiteningKeyLSB & 0x01 ) ^ ( ( WhiteningKeyLSB >> 5 ) & 0x01 );
 WhiteningKeyLSB= ( ( WhiteningKeyLSB >> 1 ) & 0xFF ) | ( ( WhiteningKeyMSBPrevious << 7 ) & 0x80 );
 }
 for( j = 1; j < bufferSize; j++ )
 {
 buffer[j] ^= WhiteningKeyLSB;

 for( i = 0; i < 8; i++ )
 {
 WhiteningKeyMSBPrevious = WhiteningKeyMSB;
 WhiteningKeyMSB = ( WhiteningKeyLSB & 0x01 ) ^ ( ( WhiteningKeyLSB >> 5 ) & 0x01 );
 WhiteningKeyLSB= ( ( WhiteningKeyLSB >> 1 ) & 0xFF ) | ( ( WhiteningKeyMSBPrevious << 7 ) & 0x80 );
 }
 }
}
*/
/*
void SX1232RadioComputeWhitening( uint8_t *buffer, uint16_t bufferSize )
{
 uint8_t i = 0;
 uint16_t j = 0;
 uint8_t WhiteningKeyMSBPrevious = 0; // 9th bit of the LFSR

 for( j = 0; j < bufferSize; j++ ) // byte counter
 {
 buffer[j] ^= WhiteningKeyLSB; // XOR between the data and the whitening key

 for( i = 0; i < 8; i++ ) // 8-bit shift between each byte
 {
 WhiteningKeyMSBPrevious = WhiteningKeyMSB;
 WhiteningKeyMSB = ( WhiteningKeyLSB & 0x01 ) ^ ( ( WhiteningKeyLSB >> 5 ) & 0x01 );
 WhiteningKeyLSB= ( ( WhiteningKeyLSB >> 1 ) & 0xFF ) | ( ( WhiteningKeyMSBPrevious << 7 ) & 0x80 );
 }
 }
}
*/

/* 
  unsigned long time;
  time = micros();
  Serial.print("Start send: ");
  Serial.println(time);
  rf95.send(data, sizeof(data));

  for(uint16_t i =0; i<1000; i++){
      digitalWrite(AMPin, LOW);
      delayMicroseconds(50);
      digitalWrite(AMPin, HIGH);
      delayMicroseconds(50);      
  }
 
  time = micros();
  Serial.print("Finish send: ");
  Serial.println(time);
  
  rf95.waitPacketSent();

  time = micros();
  Serial.print("Finish waitPacketSent: ");
  Serial.println(time);
 */

  /*
  //Different Combination for distance and speed examples: 
  Example 1: Bw = 125 kHz, Cr = 4/5, Sf = 128chips/symbol, CRC on. Default medium range
    rf95.setSignalBandwidth(125000);
    rf95.setCodingRate4(5);
    rf95.setSpreadingFactor(7);
  Example 2: Bw = 500 kHz, Cr = 4/5, Sf = 128chips/symbol, CRC on. Fast+short range
    rf95.setSignalBandwidth(500000);
    rf95.setCodingRate4(5);
    rf95.setSpreadingFactor(7);
  Example 3: Bw = 31.25 kHz, Cr = 4/8, Sf = 512chips/symbol, CRC on. Slow+long range
    rf95.setSignalBandwidth(31250);
    rf95.setCodingRate4(8);
    rf95.setSpreadingFactor(9);
  Example 4: Bw = 125 kHz, Cr = 4/8, Sf = 4096chips/symbol, CRC on. Slow+long range
    rf95.setSignalBandwidth(125000);
    rf95.setCodingRate4(8);
    rf95.setSpreadingFactor(12); 
  */
