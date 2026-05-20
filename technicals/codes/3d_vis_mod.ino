#include <WiFi.h>
#include <AsyncTCP.h>
#include <ESPAsyncWebServer.h>
#include "webnew_fin.h"

const char* ssid = "CHAOS101";
const char* password = "chaos101";

AsyncWebServer server(80);
AsyncWebSocket ws("/ws");

const int numPoints = 180;
uint16_t visualBuffer[numPoints * 3];

byte currentByte = 0;
int bitCount = 0;

void processTRNG(int valX) {
    static int lastBit = -1;
    int bit = (valX & 1);
    if (lastBit == -1) { lastBit = bit; return; }
    if (bit != lastBit) {
        currentByte = (currentByte << 1) | lastBit;
        if (++bitCount == 8) {
            ws.textAll(String(currentByte));
            bitCount = 0;
            currentByte = 0;
        }
    }
    lastBit = -1;
}

void setup() {
    Serial.begin(115200);
    analogReadResolution(12);

    WiFi.softAP(ssid, password);
    delay(100);

    Serial.print("IP: ");
    Serial.println(WiFi.softAPIP());

    server.on("/", HTTP_GET, [](AsyncWebServerRequest* request) {
        request->send(200, "text/html", webnew_fin);
    });

    server.addHandler(&ws);
    server.begin();
}

void loop() {
    // Nessun availableForWriteAll() - era quello che bloccava su Windows/Android
    if (ws.count() > 0) {
        for (int i = 0; i < numPoints; i++) {
            uint16_t x = analogRead(A0);
            uint16_t y = analogRead(A1);
            uint16_t z = analogRead(A2);

            visualBuffer[i * 3]     = x;
            visualBuffer[i * 3 + 1] = y;
            visualBuffer[i * 3 + 2] = z;

            if (i % 5 == 0) processTRNG(x);
        }
        ws.binaryAll((uint8_t*)visualBuffer, sizeof(visualBuffer));
    }

    ws.cleanupClients();
    delay(10);
}