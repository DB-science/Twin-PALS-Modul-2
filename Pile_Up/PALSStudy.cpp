/****************************************************************************
**
**  Copyright (C) 2023 Dr. Danny Petschke and Dominik Boras
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see http://www.gnu.org/licenses/.
**
*****************************************************************************
**
**  @authors: Danny Petschke, Dominik Boras
**  @contact: danny.petschke@uni-wuerzburg.de, dominik.boras@uni-wuerzburg.de
**
*****************************************************************************/

#include <iostream>
#include <algorithm>
#include <string>
#include <math.h>
#include <random>

#include "dltpulsegenerator.h"
#include "DPulseStreamAPI.h"
#include "DetEventInfo.h"


using namespace DLifeTime;


//gobal variables initialzation:

int eventIndex = 0;

int Window200nsCounter = 0;
int TrueCoincidenceCounter = 0;
int FalseCoincidenceCounter = 0;
int FalseCoincidenceCounterBackground = 0;
int PileUpCounter = 0;
int TruePileUpCounter = 0;
int FalsePileUpCounter = 0;
int BackscatterPileUpCounter = 0;
int Backscatter1275PileUpCounter = 0;
int Double511Detection = 0;

int Classification_ID = 0;
double time_fingerprint = 0.0;
double FactorTo500mv = 2.0;

//chose your PHS windows to get a guess of true coincidence events in the resulted stream


double window_511_low = (70 / FactorTo500mv);
double window_511_high = (1061 / FactorTo500mv);

double window_1275_low = (1062 / FactorTo500mv);
double window_1275_high = (3750 / FactorTo500mv);


bool Window_511_func(double NumberOfCounts) {
	if ((NumberOfCounts >= window_511_low) && (NumberOfCounts <= window_511_high)) {
		return true;
	}
	else {
		return false;
	}
}

bool Window_1275_func(double NumberOfCounts) {
	if ((NumberOfCounts >= window_1275_low) && (NumberOfCounts <= window_1275_high)) {
		return true;
	}
	else {
		return false;
	}
}

bool lookForPileUp(double time1, double time2, double pulseWidth) {
	if ((time2 - time1) >= pulseWidth) {
		return false;
	}
	else {
		return true;
	}


}

bool lookForClosePileUp(double time1, double time2, double pulseWidth) {
	if ((time2 - time1) >= pulseWidth) {
		return false;
	}
	else {
		return true;
	}


}

#define kNumberOfBins 1024
int main() {

	int SourcePower = 25; //in uCi
	int sweep = 200; //in ns
	double PulseWidth = 15.0; //in ns
	double risingEdge = 5.3; //in ns
	// (0) defines ...
	//Output Stream for drs4



	const std::string streamFileName = "G:/StreamsDissertation/RealesSetup/Geometrie/PileUp_Pulslength15ns_Zylinder1_PMMA_Sweep" + std::to_string(sweep) + "ns_" + std::to_string(SourcePower) + "uCi.drs4DataStream";


	fstream myfile("G:/StreamsDissertation/RealesSetup/Geometrie/PileUp_Pulslength15ns_Zylinder1_PMMA_AnalyseParameter" + std::to_string(SourcePower) + "uCi.txt", ios::out);


	//Input Streams

	ifstream ifile("C:/Users/Simulationsrechner/Desktop/Geant4MeetsPulsGenerator/Streams/1275keV_Zylinder1_PMMA_probe_1.5mm_Abstand_3mm.stream", ios::binary);
	ifstream ifile511("C:/Users/Simulationsrechner/Desktop/Geant4MeetsPulsGenerator/Streams/511keV_Zylinder1_PMMA_probe_1.5mm_Abstand_3mm.stream", ios::binary);

	// Initialisiere Zufallszahlengenerator
	std::random_device rd;
	std::mt19937 generator(rd());

	// Initialisiere exponentielle Verteilungen
	std::exponential_distribution<double> distribution1(1.0);
	std::exponential_distribution<double> distribution2(1.0);
	std::exponential_distribution<double> distribution3(1.0);
	std::exponential_distribution<double> distribution4(1.0);
	std::exponential_distribution<double> distribution5(1.0);
	
	
	//positron lifetime!
	int NumberOfComponentes = 3; //the number of Components the Positron lifetime Spectrum should have. 1 to 5 Components availible, Intensity must add to 1.
	double lifetime1 = 0.158;double intensity1 = 0.825;
	double lifetime2 = 0.380;double intensity2 = 0.172;
	double lifetime3 = 2.750;double intensity3 = 0.003;
	double lifetime4 = 0.0; double intensity4 = 0.00;
	double lifetime5 = 0.0; double intensity5 = 0.00;

	

	//definition of vectors to manage events
	std::vector<int> sequence;
	std::vector<double> Det1;
	std::vector<double> Det2;
	std::vector<int> Det1numberOfCounts;
	std::vector<int> Det2numberOfCounts;
	std::vector<int> Det1event;
	std::vector<int> Det2event;

	std::vector<int> Det1eventValid;
	std::vector<int> Det2eventValid;

	std::vector<bool> Det1valid;
	std::vector<bool> Det2valid;

	std::vector<std::pair<double, double> > A, B; // example code!



	double activity = 37000 * SourcePower;

	// Calculate the decay constant based on the activity level
	double decay_constant = 1 / activity;

	// Create a random number generator
	std::random_device rd1;
	std::mt19937 gen2(rd1());

	std::exponential_distribution<double> distribution(1 / decay_constant);


	//DLT -Setup:

	DLTSetup setup = DLTSetup_DEMO;

	setup.sweep = sweep; // ns
	setup.numberOfCells = kNumberOfBins; // 1024 # capacitors
	setup.ATS = 0.; // arrival time spread (you can leave it as it is)

	DLTPulse pulse = DLTPulse2_DEMO;

	pulse.amplitude = -500.; // mV --> scale all amplitudes into this range!
	pulse.isPositiveSignalPolarity = false; // negative pulses
	pulse.delay = 0.; // trigger delay (you can leave it as it is)

	// detector A:

	pulse.pulseA.riseTime = 5.3;   // ns
	pulse.pulseA.pulseWidth = 0.315; // ns
	pulse.pulseA.randomNoiseInfoV.enabled = true;
	pulse.pulseA.randomNoiseInfoV.rndNoise = 0.35;  // rms [mW]

	pulse.pulseA.baselineOffsetJitterInfoV.enabled = true;
	pulse.pulseA.baselineOffsetJitterInfoV.meanOfBaselineOffsetJitter = 0.; // mV
	pulse.pulseA.baselineOffsetJitterInfoV.stddevOfBaselineOffsetJitter = 5.; // mV

	pulse.pulseA.timeAxisNonLinearityInfoT.enabled = true;
	pulse.pulseA.timeAxisNonLinearityInfoT.fixedPatternApertureJitter = 0.005;  // ns
	pulse.pulseA.timeAxisNonLinearityInfoT.rndApertureJitter = 0.0025; // ns

	// detector B:

	pulse.pulseB.riseTime = 5.3;   // ns
	pulse.pulseB.pulseWidth = 0.315; // ns
	pulse.pulseB.randomNoiseInfoV.enabled = true;
	pulse.pulseB.randomNoiseInfoV.rndNoise = 0.35;  // rms [mW]

	pulse.pulseB.baselineOffsetJitterInfoV.enabled = true;
	pulse.pulseB.baselineOffsetJitterInfoV.meanOfBaselineOffsetJitter = 0.;   // mV
	pulse.pulseB.baselineOffsetJitterInfoV.stddevOfBaselineOffsetJitter = 5.; // mV

	pulse.pulseB.timeAxisNonLinearityInfoT.enabled = true;
	pulse.pulseB.timeAxisNonLinearityInfoT.fixedPatternApertureJitter = 0.005;  // ns
	pulse.pulseB.timeAxisNonLinearityInfoT.rndApertureJitter = 0.0025;          // ns

	// unused for this purpose but it needs to be defined ...

	DLTPHS phs = DLTPHS2_DEMO;            // not used! 
	DLTSimulationInput simulationInput = DLTSimulationInput_DEMO; // not used!

	phs.useGaussianModels = true;
	// (A) generate pulsegenerator instance ...

	CallbackClass callback; // check for errors ...

	DLTPulseGenerator* pulseGenerator = new DLTPulseGenerator(simulationInput, phs, setup, pulse, &callback);

	if (callback.hasError()) {
		DDELETE_SAFETY(pulseGenerator);

		std::cout << "error: pulse generator cannot be initialized successfully!";

		return 0;
	}

	// (B) create pulse stream file ...

	if (!DPulseStreamManager::sharedInstance()->start(streamFileName, setup.sweep, double(kNumberOfBins) / setup.sweep, kNumberOfBins)) { // DRS4 compatible ...
		std::cout << "file cannot be created!";

		return 0;
	}

	const int sizeOfWave = sizeof(float) * kNumberOfBins;


	for (int i = 0; i < 50; i++) {
		int stream2Index = i + 1; // Leseindex für stream2
		int stream1Index = 0;


		if (ifile.is_open() || ifile511.is_open())
		{


			while (!ifile.eof()) {

				if (stream2Index >= 50000000) {
					stream2Index = 0;
				}

				EventInfo info;
				EventInfo info511;

				DLTPulseF pulseA, pulseB;


				ifile.seekg(stream1Index * sizeof(EventInfo));
				ifile.read(reinterpret_cast<char*>(&info), sizeof(EventInfo));
				++stream1Index;

				ifile.read((char*)&info, sizeof(EventInfo));

				ifile511.seekg(stream2Index * sizeof(EventInfo));
				ifile511.read(reinterpret_cast<char*>(&info511), sizeof(EventInfo));
				++stream2Index;




				//definition of geant4stream variables:
				double NumberOfCounts1 = (info.info1().numberOfCounts() / FactorTo500mv);
				double NumberOfCounts2 = (info.info2().numberOfCounts() / FactorTo500mv);

				double NumberOfCounts1_511 = (info511.info1().numberOfCounts() / FactorTo500mv);
				double NumberOfCounts2_511 = (info511.info2().numberOfCounts() / FactorTo500mv);

				double GlTime1274 = info.startTime();
				double GlTime511 = info511.startTime();

				double PALS = positronAnhillationLifetime(NumberOfComponentes,lifetime1, lifetime2, lifetime3, lifetime4, lifetime5,
														intensity1, intensity2, intensity3, intensity4, intensity5,
														distribution1, distribution2, distribution3, distribution4, distribution5,
														generator);


				double pmt1 = generatePMTuncertaincy(0.55, 1.0, info.info1().numberOfCounts());

				double pmt1_511 = generatePMTuncertaincy(0.55, 1.0, info511.info1().numberOfCounts());

				double pmt2 = generatePMTuncertaincy(0.55, 1.0, info.info2().numberOfCounts());

				double pmt2_511 = generatePMTuncertaincy(0.55, 1.0, info511.info2().numberOfCounts());

				//1275 Gamma time
				double Det1_1275_arrivalTime = 5 + (info.info1().MeanArr() - GlTime1274) + pmt1;

				double Det2_1275_arrivalTime = 5 + (info.info2().MeanArr() - GlTime1274) + pmt2;

				//511 Gamma time
				double Det1_511_arrivalTime = 5 + (info511.info1().MeanArr() - GlTime511) + pmt1_511 + PALS;//

				double Det2_511_arrivalTime = 5 + (info511.info2().MeanArr() - GlTime511) + pmt2_511 + PALS;//



				int event_ID = info.info1().id();


				//classification true coincidence:

				const bool Coincidence1 = ((info.info1().isValid() && info511.info2().isValid()) && ((info.info2().isValid() == false) && (info511.info1().isValid() == false)));
				const bool Coincidence2 = ((info.info2().isValid() && info511.info1().isValid()) && ((info.info1().isValid() == false) && (info511.info2().isValid() == false)));

				//classification bad events:
				const bool Bad511 = (info511.info1().isValid() && info511.info2().isValid()) && (info.info1().isValid() == false && info.info2().isValid() == false);
				const bool Bad1275 = (info.info1().isValid() && info.info2().isValid()) && ((info511.info1().isValid() == false && info511.info2().isValid() == false));

				//Backscattering Coincidence
				const bool Bad511Coincidence1 = ((info511.info1().isValid() && info511.info2().isValid()) && (info.info1().isValid() && info.info2().isValid() == false));
				const bool Bad511Coincidence2 = ((info511.info1().isValid() && info511.info2().isValid()) && (info.info2().isValid() && info.info1().isValid() == false));

				const bool Bad1275Coincidence1 = ((info.info1().isValid() && info.info2().isValid()) && (info511.info1().isValid() && info511.info2().isValid() == false));
				const bool Bad1275Coincidence2 = ((info.info1().isValid() && info.info2().isValid()) && (info511.info2().isValid() && info511.info1().isValid() == false));

				const bool WorstCase = (info.info1().isValid() && info.info2().isValid()) && (info511.info1().isValid() && info511.info2().isValid());


				//single events

				const bool single1_511 = (info511.info1().isValid() && ((info511.info2().isValid() == false) && (info.info1().isValid() == false) && (info.info2().isValid() == false)));
				const bool single2_511 = (info511.info2().isValid() && ((info511.info1().isValid() == false) && (info.info1().isValid() == false) && (info.info2().isValid() == false)));

				const bool single1_1275 = (info.info1().isValid() && ((info.info2().isValid() == false) && (info511.info1().isValid() == false) && (info511.info2().isValid() == false)));
				const bool single2_1275 = (info.info2().isValid() && ((info.info1().isValid() == false) && (info511.info1().isValid() == false) && (info511.info2().isValid() == false)));

				double Source_decay = radioactivity(decay_constant, distribution, gen2);

				eventIndex++;

				// generate pulses from classified Geant4 events ...
				if (time_fingerprint >= sweep) {
					if (((A.size() >= 1 && B.size() > 1) && ((Det2[Det2.size() - 1] - Det2[Det2.size() - 2]) < 15)) || ((A.size() > 1 && B.size() >= 1) && ((Det1[Det1.size() - 1] - Det1[Det1.size() - 2]) < 15)))
					{

						if ((A.size() == 1) && (B.size() == 1)) {
							if (Det1eventValid[0] == Det2eventValid[0]) {
								if (Det1valid[0] == true && Det2valid[0] == true) {
									if (sequence[0] == 1 || sequence[0] == 2) {
										if ((Window_1275_func(Det1numberOfCounts[0]) && Window_511_func(Det2numberOfCounts[0])) || (Window_1275_func(Det2numberOfCounts[0]) && Window_511_func(Det1numberOfCounts[0])))
										{
											Window200nsCounter++;
											TrueCoincidenceCounter++;
										}
									}

									if (sequence[0] == 7) {
										if ((Window_1275_func(Det1numberOfCounts[0]) && Window_511_func(Det2numberOfCounts[0])) || (Window_1275_func(Det2numberOfCounts[0]) && Window_511_func(Det1numberOfCounts[0])))
										{
											Window200nsCounter++;
											Backscatter1275PileUpCounter++;
										}
									}

									if (sequence[0] == 8) {
										Window200nsCounter++;
										Double511Detection++;
									}

								}
							}
							if (Det1eventValid[0] != Det2eventValid[0]) {
								if (Det1valid[0] == true && Det2valid[0] == true) {
									if ((Window_1275_func(Det1numberOfCounts[0]) && Window_511_func(Det2numberOfCounts[0])))
									{
										double time = Det2[0] - Det1[0];
										if (-8.0 <= time <= 42.0) {
											Window200nsCounter++;
											FalseCoincidenceCounter++;
										}

									}
									if ((Window_1275_func(Det2numberOfCounts[0]) && Window_511_func(Det1numberOfCounts[0]))) {
										double time = Det1[0] - Det2[0];
										if (-8.0 <= time <= 42.0) {
											Window200nsCounter++;
											FalseCoincidenceCounter++;
										}
									}
								}
							}


						}

						if ((A.size() + B.size()) >= 3) {
							if (Det1eventValid[0] == Det2eventValid[0]) {
								if ((Window_1275_func(Det1numberOfCounts[0]) && Window_511_func(Det2numberOfCounts[0])) || (Window_1275_func(Det2numberOfCounts[0]) && Window_511_func(Det1numberOfCounts[0])))
								{
									Window200nsCounter++;
									TruePileUpCounter++;
								}
							}

							if (Det1eventValid[0] != Det2eventValid[0]) {
								if ((Window_1275_func(Det1numberOfCounts[0]) && Window_511_func(Det2numberOfCounts[0])))
								{
									double time = Det2[0] - Det1[0];
									if (-8.0 <= time <= 42.0) {
										Window200nsCounter++;
										FalsePileUpCounter++;
									}

								}
								if ((Window_1275_func(Det2numberOfCounts[0]) && Window_511_func(Det1numberOfCounts[0]))) {
									double time = Det1[0] - Det2[0];
									if (-8.0 <= time <= 42.0) {
										Window200nsCounter++;
										FalsePileUpCounter++;
									}
								}

							}
						}

						if (pulseGenerator->emitPulsesByExternalInput(&pulseA, &pulseB, A, B)) {
							float tChannel0[kNumberOfBins] = { 0 };
							float tChannel1[kNumberOfBins] = { 0 };

							float waveChannel0[kNumberOfBins] = { 0 };
							float waveChannel1[kNumberOfBins] = { 0 };

							std::fill(tChannel0, tChannel0 + kNumberOfBins, 0.);
							std::fill(tChannel1, tChannel1 + kNumberOfBins, 0.);

							std::fill(waveChannel0, waveChannel0 + kNumberOfBins, 0.);
							std::fill(waveChannel1, waveChannel1 + kNumberOfBins, 0.);

							for (int i = 0; i < kNumberOfBins; ++i) {
								tChannel0[i] = pulseA.at(i).x();
								waveChannel0[i] = pulseA.at(i).y();

								tChannel1[i] = pulseB.at(i).x();
								waveChannel1[i] = pulseB.at(i).y();
							}
							// append pulses to stream ...
							DPulseStreamManager::sharedInstance()->writePulsePair(tChannel0, waveChannel0, tChannel1, waveChannel1, sizeOfWave);

						}

						A.clear();
						B.clear();

						time_fingerprint = 0.0;
						Det1numberOfCounts.clear();
						Det2numberOfCounts.clear();
						Det1.clear();
						Det2.clear();
						Det1event.clear();
						Det2event.clear();
						sequence.clear();
						Det1valid.clear();
						Det2valid.clear();
						Det1eventValid.clear();
						Det2eventValid.clear();
					}
					A.clear();
					B.clear();

					time_fingerprint = 0.0;
					Det1numberOfCounts.clear();
					Det2numberOfCounts.clear();
					Det1.clear();
					Det2.clear();
					Det1event.clear();
					Det2event.clear();
					sequence.clear();
					Det1valid.clear();
					Det2valid.clear();
					Det1eventValid.clear();
					Det2eventValid.clear();
				}

				
				if (Coincidence1) {
					Classification_ID = 1;
					sequence.push_back(Classification_ID);
					Det1.push_back(Det1_1275_arrivalTime + time_fingerprint);
					Det2.push_back(Det2_511_arrivalTime + time_fingerprint);
					Det1numberOfCounts.push_back(NumberOfCounts1);
					Det2numberOfCounts.push_back(NumberOfCounts2_511);
					Det1event.push_back(eventIndex);
					Det2event.push_back(eventIndex);

					if (Window_1275_func(NumberOfCounts1) || Window_511_func(NumberOfCounts1)) {
						Det1valid.push_back(true);
					}
					else {
						Det1valid.push_back(false);
					}

					if (Window_511_func(NumberOfCounts2_511) || Window_511_func(NumberOfCounts2_511)) {
						Det2valid.push_back(true);
					}
					else {
						Det2valid.push_back(false);
					}
					if (Det1.size() > 1) {
						if ((Det1valid[Det1valid.size() - 1] == true) && (Det1valid[Det1valid.size() - 2] == false) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det1[Det1.size() - 2]) , -(Det1numberOfCounts[Det1numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back((Det1event[Det1event.size() - 2]));

							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

						if ((Det1valid[Det1valid.size() - 2] == true) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_1275_func(NumberOfCounts1) || Window_511_func(NumberOfCounts1)) {
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}
					}
					if (Det2.size() > 1) {

						if ((Det2valid[Det2valid.size() - 1] == true) && (Det2valid[Det2valid.size() - 2] == false) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det2[Det2.size() - 2]) , -(Det2numberOfCounts[Det2numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det2eventValid.push_back((Det2event[Det2event.size() - 2]));
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

						if ((Det2valid[Det2valid.size() - 2] == true) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_511_func(NumberOfCounts2_511) || Window_1275_func(NumberOfCounts2_511)) {
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}
					}
				}

				if (Coincidence2) {
					Classification_ID = 2;
					sequence.push_back(Classification_ID);
					Det1.push_back(Det1_511_arrivalTime + time_fingerprint);
					Det2.push_back(Det2_1275_arrivalTime + time_fingerprint);
					Det1numberOfCounts.push_back(NumberOfCounts1_511);
					Det2numberOfCounts.push_back(NumberOfCounts2);
					Det1event.push_back(eventIndex);
					Det2event.push_back(eventIndex);

					if (Window_1275_func(NumberOfCounts1_511) || Window_511_func(NumberOfCounts1_511)) {
						Det1valid.push_back(true);
					}
					else {
						Det1valid.push_back(false);
					}

					if (Window_511_func(NumberOfCounts2) || Window_1275_func(NumberOfCounts2)) {
						Det2valid.push_back(true);
					}
					else {
						Det2valid.push_back(false);
					}
					if (Det1.size() > 1) {

						if ((Det1valid[Det1valid.size() - 1] == true) && (Det1valid[Det1valid.size() - 2] == false) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det1[Det1.size() - 2]) , -(Det1numberOfCounts[Det1numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back((Det1event[Det1event.size() - 2]));
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

						if ((Det1valid[Det1valid.size() - 2] == true) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}


					}
					else {
						if (Window_1275_func(NumberOfCounts1_511) || Window_511_func(NumberOfCounts1_511)) {
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}
					}
					if (Det2.size() > 1) {
						if ((Det2valid[Det2valid.size() - 1] == true) && (Det2valid[Det2valid.size() - 2] == false) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det2[Det2.size() - 2]) , -(Det2numberOfCounts[Det2numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det2eventValid.push_back((Det2event[Det2event.size() - 2]));
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}
						if ((Det2valid[Det2valid.size() - 2] == true) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_511_func(NumberOfCounts2) || Window_1275_func(NumberOfCounts2)) {
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}
					}

				}
				
				if (single1_511) {
					Classification_ID = 3;
					sequence.push_back(Classification_ID);
					Det1.push_back(Det1_511_arrivalTime + time_fingerprint);
					Det1numberOfCounts.push_back(NumberOfCounts1_511);
					Det1event.push_back(eventIndex);

					if (Window_1275_func(NumberOfCounts1_511) || Window_511_func(NumberOfCounts1_511)) {
						Det1valid.push_back(true);
					}
					else {
						Det1valid.push_back(false);
					}

					if (Det1.size() > 1) {

						if ((Det1valid[Det1valid.size() - 1] == true) && (Det1valid[Det1valid.size() - 2] == false) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det1[Det1.size() - 2]) , -(Det1numberOfCounts[Det1numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back((Det1event[Det1event.size() - 2]));
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

						if ((Det1valid[Det1valid.size() - 2] == true) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}


					}
					else {
						if (Window_1275_func(NumberOfCounts1_511) || Window_511_func(NumberOfCounts1_511)) {
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}
					}
				}

				if (single2_511) {
					Classification_ID = 4;
					sequence.push_back(Classification_ID);
					Det2.push_back(Det2_511_arrivalTime + time_fingerprint);
					Det2numberOfCounts.push_back(NumberOfCounts2_511);
					Det2event.push_back(eventIndex);

					if (Window_511_func(NumberOfCounts2_511) || Window_1275_func(NumberOfCounts2_511)) {
						Det2valid.push_back(true);
					}
					else {
						Det2valid.push_back(false);
					}

					if (Det2.size() > 1) {

						if ((Det2valid[Det2valid.size() - 1] == true) && (Det2valid[Det2valid.size() - 2] == false) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det2[Det2.size() - 2]) , -(Det2numberOfCounts[Det2numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det2eventValid.push_back((Det2event[Det2event.size() - 2]));
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

						if ((Det2valid[Det2valid.size() - 2] == true) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_511_func(NumberOfCounts2_511) || Window_1275_func(NumberOfCounts2_511)) {
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}
					}
				}

				if (single1_1275) {
					Classification_ID = 5;
					sequence.push_back(Classification_ID);
					Det1.push_back(Det1_1275_arrivalTime + time_fingerprint);
					Det1numberOfCounts.push_back(NumberOfCounts1);
					Det1event.push_back(eventIndex);

					if (Window_1275_func(NumberOfCounts1) || Window_511_func(NumberOfCounts1)) {
						Det1valid.push_back(true);
					}
					else {
						Det1valid.push_back(false);
					}

					if (Det1.size() > 1) {
						if ((Det1valid[Det1valid.size() - 1] == true) && (Det1valid[Det1valid.size() - 2] == false) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det1[Det1.size() - 2]) , -(Det1numberOfCounts[Det1numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back((Det1event[Det1event.size() - 2]));
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

						if ((Det1valid[Det1valid.size() - 2] == true) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_1275_func(NumberOfCounts1) || Window_511_func(NumberOfCounts1)) {
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}
					}

				}
				if (single2_1275) {
					Classification_ID = 6;
					sequence.push_back(Classification_ID);
					Det2.push_back(Det2_1275_arrivalTime + time_fingerprint);
					Det2numberOfCounts.push_back(NumberOfCounts2);
					Det2event.push_back(eventIndex);

					if (Window_511_func(NumberOfCounts2) || Window_1275_func(NumberOfCounts2)) {
						Det2valid.push_back(true);
					}
					else {
						Det2valid.push_back(false);
					}

					if (Det2.size() > 1) {
						if ((Det2valid[Det2valid.size() - 1] == true) && (Det2valid[Det2valid.size() - 2] == false) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det2[Det2.size() - 2]) , -(Det2numberOfCounts[Det2numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det2eventValid.push_back((Det2event[Det2event.size() - 2]));
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}
						if ((Det2valid[Det2valid.size() - 2] == true) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_511_func(NumberOfCounts2) || Window_1275_func(NumberOfCounts2)) {
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}
					}
				}
				
				if (Bad1275) {
					Classification_ID = 7;
					sequence.push_back(Classification_ID);
					Det1.push_back(Det1_1275_arrivalTime + time_fingerprint);
					Det2.push_back(Det2_1275_arrivalTime + time_fingerprint);
					Det1event.push_back(eventIndex);
					Det2event.push_back(eventIndex);
					Det1numberOfCounts.push_back(NumberOfCounts1);
					Det2numberOfCounts.push_back(NumberOfCounts2);

					if (Window_1275_func(NumberOfCounts1) || Window_511_func(NumberOfCounts1)) {
						Det1valid.push_back(true);
					}
					else {
						Det1valid.push_back(false);
					}

					if (Window_511_func(NumberOfCounts2) || Window_1275_func(NumberOfCounts2)) {
						Det2valid.push_back(true);
					}
					else {
						Det2valid.push_back(false);
					}

					if (Det1.size() > 1) {

						if ((Det1valid[Det1valid.size() - 1] == true) && (Det1valid[Det1valid.size() - 2] == false) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det1[Det1.size() - 2]) , -(Det1numberOfCounts[Det1numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back((Det1event[Det1event.size() - 2]));
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

						if ((Det1valid[Det1valid.size() - 2] == true) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_1275_func(NumberOfCounts1) || Window_511_func(NumberOfCounts1)) {
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}
					}

					if (Det2.size() > 1) {

						if ((Det2valid[Det2valid.size() - 1] == true) && (Det2valid[Det2valid.size() - 2] == false) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det2[Det2.size() - 2]) , -(Det2numberOfCounts[Det2numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det2eventValid.push_back((Det2event[Det2event.size() - 2]));
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}
						if ((Det2valid[Det2valid.size() - 2] == true) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_511_func(NumberOfCounts2) || Window_1275_func(NumberOfCounts2)) {
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}
					}

				}
				
				if (Bad511) {
					Classification_ID = 8;
					sequence.push_back(Classification_ID);
					Det1.push_back(Det1_511_arrivalTime + time_fingerprint);
					Det2.push_back(Det2_511_arrivalTime + time_fingerprint);
					Det1event.push_back(eventIndex);
					Det2event.push_back(eventIndex);
					Det1numberOfCounts.push_back(NumberOfCounts1_511);
					Det2numberOfCounts.push_back(NumberOfCounts2_511);

					if (Window_1275_func(NumberOfCounts1_511) || Window_511_func(NumberOfCounts1_511)) {
						Det1valid.push_back(true);
					}
					else {
						Det1valid.push_back(false);
					}

					if (Window_511_func(NumberOfCounts2_511) || Window_1275_func(NumberOfCounts2_511)) {
						Det2valid.push_back(true);
					}
					else {
						Det2valid.push_back(false);
					}

					if (Det1.size() > 1) {

						if ((Det1valid[Det1valid.size() - 1] == true) && (Det1valid[Det1valid.size() - 2] == false) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det1[Det1.size() - 2]) , -(Det1numberOfCounts[Det1numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back((Det1event[Det1event.size() - 2]));
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

						if ((Det1valid[Det1valid.size() - 2] == true) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_1275_func(NumberOfCounts1_511) || Window_511_func(NumberOfCounts1_511)) {
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}
					}

					if (Det2.size() > 1) {

						if ((Det2valid[Det2valid.size() - 1] == true) && (Det2valid[Det2valid.size() - 2] == false) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det2[Det2.size() - 2]) , -(Det2numberOfCounts[Det2numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det2eventValid.push_back((Det2event[Det2event.size() - 2]));
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

						if ((Det2valid[Det2valid.size() - 2] == true) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_511_func(NumberOfCounts2_511) || Window_1275_func(NumberOfCounts2_511)) {
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}
					}



				}

				/*
				if (Bad511Coincidence1) {
					Classification_ID = 9;
					sequence.push_back(Classification_ID);
					//Det1.push_back(Det1_511_arrivalTime+time_fingerprint);
					Det2.push_back(Det2_511_arrivalTime + time_fingerprint);
					Det1.push_back(Det1_1275_arrivalTime + time_fingerprint);
					Det1event.push_back(eventIndex);
					Det2event.push_back(eventIndex);
					Det1numberOfCounts.push_back(NumberOfCounts1);
					Det2numberOfCounts.push_back(NumberOfCounts2_511);

					if ((NumberOfCounts1 + (0.5 * NumberOfCounts1_511)) >= (window_511_low)) {
						Det1valid.push_back(true);
					}
					else {
						Det1valid.push_back(false);
					}

					if (Window_511_func(NumberOfCounts2_511) || Window_1275_func(NumberOfCounts2_511)) {
						Det2valid.push_back(true);
					}
					else {
						Det2valid.push_back(false);
					}

					if (Det1.size() > 1) {

						if ((Det1valid[Det1valid.size() - 1] == true) && (Det1valid[Det1valid.size() - 2] == false) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det1[Det1.size() - 2]) , -(Det1numberOfCounts[Det1numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back((Det1event[Det1event.size() - 2]));
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(d);
							Det1eventValid.push_back(eventIndex);
							Det1.push_back(Det1_511_arrivalTime + time_fingerprint);
							Det1numberOfCounts.push_back(NumberOfCounts1_511);
						}

						if ((Det1valid[Det1valid.size() - 2] == true) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(d);
							Det1eventValid.push_back(eventIndex);
							Det1.push_back(Det1_511_arrivalTime + time_fingerprint);
							Det1numberOfCounts.push_back(NumberOfCounts1_511);
						}

					}
					else {
						if ((NumberOfCounts1 + (0.5 * NumberOfCounts1_511)) >= (window_511_low)) {
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(d);
							Det1eventValid.push_back(eventIndex);
							Det1.push_back(Det1_511_arrivalTime + time_fingerprint);
							Det1numberOfCounts.push_back(NumberOfCounts1_511);
						}
					}

					if (Det2.size() > 1) {

						if ((Det2valid[Det2valid.size() - 1] == true) && (Det2valid[Det2valid.size() - 2] == false) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det2[Det2.size() - 2]) , -(Det2numberOfCounts[Det2numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							B.push_back(b);
							Det2eventValid.push_back((Det2event[Det2event.size() - 2]));
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

						if ((Det2valid[Det2valid.size() - 2] == true) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_511_func(NumberOfCounts2_511) || Window_1275_func(NumberOfCounts2_511)) {
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}
					}
				}


				if (Bad511Coincidence2) {
					Classification_ID = 10;
					sequence.push_back(Classification_ID);
					Det1.push_back(Det1_511_arrivalTime + time_fingerprint);
					//Det2.push_back(IRF_511_2+time_fingerprint);
					Det2.push_back(Det2_1275_arrivalTime + time_fingerprint);
					Det1event.push_back(eventIndex);
					Det2event.push_back(eventIndex);
					Det1numberOfCounts.push_back(NumberOfCounts1_511);
					Det2numberOfCounts.push_back(NumberOfCounts2);

					if (Window_1275_func(NumberOfCounts1_511) || Window_511_func(NumberOfCounts1_511)) {
						Det1valid.push_back(true);
					}
					else {
						Det1valid.push_back(false);
					}

					if ((NumberOfCounts2 + (0.5 * NumberOfCounts2_511)) >= (window_511_low)) {
						Det2valid.push_back(true);
					}
					else {
						Det2valid.push_back(false);
					}

					if (Det1.size() > 1) {

						if ((Det1valid[Det1valid.size() - 1] == true) && (Det1valid[Det1valid.size() - 2] == false) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det1[Det1.size() - 2]) , -(Det1numberOfCounts[Det1numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back((Det1event[Det1event.size() - 2]));
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

						if ((Det1valid[Det1valid.size() - 2] == true) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_1275_func(NumberOfCounts1_511) || Window_511_func(NumberOfCounts1_511)) {
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}
					}

					if (Det2.size() > 1) {

						if ((Det2valid[Det2valid.size() - 1] == true) && (Det2valid[Det2valid.size() - 2] == false) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det2[Det2.size() - 2]) , -(Det2numberOfCounts[Det2numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							B.push_back(b);
							Det2eventValid.push_back((Det2event[Det2event.size() - 2]));
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(d);
							Det2eventValid.push_back(eventIndex);
							Det2.push_back(Det2_511_arrivalTime + time_fingerprint);
							Det2numberOfCounts.push_back(NumberOfCounts2_511);
						}

						if ((Det2valid[Det2valid.size() - 2] == true) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(d);
							Det2eventValid.push_back(eventIndex);
							Det2.push_back(Det2_511_arrivalTime + time_fingerprint);
							Det2numberOfCounts.push_back(NumberOfCounts2_511);
						}

					}
					else {
						if ((NumberOfCounts2 + (0.5 * NumberOfCounts2_511)) >= (window_511_low)) {
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(d);
							Det2eventValid.push_back(eventIndex);
							Det2.push_back(Det2_511_arrivalTime + time_fingerprint);
							Det2numberOfCounts.push_back(NumberOfCounts2_511);
						}
					}


				}

				if (Bad1275Coincidence1) {
					Classification_ID = 11;
					sequence.push_back(Classification_ID);
					//Det1.push_back(IRF_1274_1+time_fingerprint);
					Det2.push_back(Det2_1275_arrivalTime + time_fingerprint);
					Det1.push_back(Det1_511_arrivalTime + time_fingerprint);
					Det1event.push_back(eventIndex);
					Det2event.push_back(eventIndex);
					Det1numberOfCounts.push_back(NumberOfCounts1_511);
					Det2numberOfCounts.push_back(NumberOfCounts2);


					if ((NumberOfCounts1 + (0.5 * NumberOfCounts1_511)) >= (window_511_low)) {
						Det1valid.push_back(true);
					}
					else {
						Det1valid.push_back(false);
					}

					if (Window_511_func(NumberOfCounts2) || Window_1275_func(NumberOfCounts2)) {
						Det2valid.push_back(true);
					}
					else {
						Det2valid.push_back(false);
					}

					if (Det1.size() > 1) {

						if ((Det1valid[Det1valid.size() - 1] == true) && (Det1valid[Det1valid.size() - 2] == false) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det1[Det1.size() - 2]) , -(Det1numberOfCounts[Det1numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back((Det1event[Det1event.size() - 2]));
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(d);
							Det1eventValid.push_back(eventIndex);
							Det1.push_back(Det1_1275_arrivalTime + time_fingerprint);
							Det1numberOfCounts.push_back(NumberOfCounts1);

						}

						if ((Det1valid[Det1valid.size() - 2] == true) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(d);
							Det1eventValid.push_back(eventIndex);
							Det1.push_back(Det1_1275_arrivalTime + time_fingerprint);
							Det1numberOfCounts.push_back(NumberOfCounts1);
						}

					}
					else {
						if ((NumberOfCounts1 + (0.5 * NumberOfCounts1_511)) >= (window_511_low)) {
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(d);
							Det1eventValid.push_back(eventIndex);
							Det1.push_back(Det1_1275_arrivalTime + time_fingerprint);
							Det1numberOfCounts.push_back(NumberOfCounts1);
						}
					}

					if (Det2.size() > 1) {

						if ((Det2valid[Det2valid.size() - 1] == true) && (Det2valid[Det2valid.size() - 2] == false) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det2[Det2.size() - 2]) , -(Det2numberOfCounts[Det2numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							B.push_back(b);
							Det2eventValid.push_back((Det2event[Det2event.size() - 2]));
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

						if ((Det2valid[Det2valid.size() - 2] == true) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_511_func(NumberOfCounts2) || Window_1275_func(NumberOfCounts2)) {
							std::pair<double, double> c = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
						}
					}

				}

				if (Bad1275Coincidence2) {
					Classification_ID = 12;
					sequence.push_back(Classification_ID);
					Det1.push_back(Det1_1275_arrivalTime + time_fingerprint);
					//Det2.push_back(IRF_1274_2+time_fingerprint);
					Det2.push_back(Det2_511_arrivalTime + time_fingerprint);
					Det1event.push_back(eventIndex);
					Det2event.push_back(eventIndex);
					Det1numberOfCounts.push_back(NumberOfCounts1);
					Det2numberOfCounts.push_back(NumberOfCounts2_511);

					if (Window_1275_func(NumberOfCounts1) || Window_511_func(NumberOfCounts1)) {
						Det1valid.push_back(true);
					}
					else {
						Det1valid.push_back(false);
					}

					if ((NumberOfCounts2 + (0.5 * NumberOfCounts2_511)) >= (window_511_low)) {
						Det2valid.push_back(true);
					}
					else {
						Det2valid.push_back(false);
					}

					if (Det1.size() > 1) {

						if ((Det1valid[Det1valid.size() - 1] == true) && (Det1valid[Det1valid.size() - 2] == false) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det1[Det1.size() - 2]) , -(Det1numberOfCounts[Det1numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back((Det1event[Det1event.size() - 2]));
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

						if ((Det1valid[Det1valid.size() - 2] == true) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}

					}
					else {
						if (Window_1275_func(NumberOfCounts1) || Window_511_func(NumberOfCounts1)) {
							std::pair<double, double> a = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
						}
					}

					if (Det2.size() > 1) {

						if ((Det2valid[Det2valid.size() - 1] == true) && (Det2valid[Det2valid.size() - 2] == false) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det2[Det2.size() - 2]) , -(Det2numberOfCounts[Det2numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							B.push_back(b);
							Det2eventValid.push_back((Det2event[Det2event.size() - 2]));
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
							Det2.push_back(Det2_1275_arrivalTime + time_fingerprint);
							Det2numberOfCounts.push_back(NumberOfCounts2);
							std::pair<double, double> d = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(d);
							Det2eventValid.push_back(eventIndex);

						}

						if ((Det2valid[Det2valid.size() - 2] == true) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(d);
							Det2eventValid.push_back(eventIndex);
							Det2.push_back(Det2_1275_arrivalTime + time_fingerprint);
							Det2numberOfCounts.push_back(NumberOfCounts2);
						}

					}
					else {
						if ((NumberOfCounts2 + (0.5 * NumberOfCounts2_511)) >= (window_511_low)) {
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(d);
							Det2eventValid.push_back(eventIndex);
							Det2.push_back(Det2_1275_arrivalTime + time_fingerprint);
							Det2numberOfCounts.push_back(NumberOfCounts2);
						}
					}

				}

				if (WorstCase) {

					Classification_ID = 13;
					sequence.push_back(Classification_ID);
					Det1.push_back(Det1_1275_arrivalTime + time_fingerprint);
					Det2.push_back(Det2_1275_arrivalTime + time_fingerprint);
					//Det2.push_back(Det2_511_arrivalTime + time_fingerprint);
					//Det1.push_back(Det1_511_arrivalTime + time_fingerprint);
					Det1event.push_back(eventIndex);
					Det2event.push_back(eventIndex);
					//Det1event.push_back(eventIndex);
					//Det2event.push_back(eventIndex);
					Det1numberOfCounts.push_back(NumberOfCounts1);
					//Det2numberOfCounts.push_back(NumberOfCounts2_511);
					//Det1numberOfCounts.push_back(NumberOfCounts1_511);
					Det2numberOfCounts.push_back(NumberOfCounts2);

					if ((NumberOfCounts1 + (0.5 * NumberOfCounts1_511)) >= (window_511_low)) {
						Det1valid.push_back(true);
					}
					else {
						Det1valid.push_back(false);

					}

					if ((NumberOfCounts2 + (0.5 * NumberOfCounts2_511)) >= (window_511_low)) {
						Det2valid.push_back(true);
					}
					else {
						Det2valid.push_back(false);
					}

					if (Det1.size() > 2) {

						if ((Det1valid[Det1valid.size() - 1] == true) && (Det1valid[Det1valid.size() - 2] == false) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {

							std::pair<double, double> b = { (Det1[Det1.size() - 2]) , -(Det1numberOfCounts[Det1numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back((Det1event[Det1event.size() - 2]));

							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(d);
							Det1eventValid.push_back(eventIndex);


						}

						if ((Det1valid[Det1valid.size() - 2] == true) && lookForPileUp(Det1[Det1.size() - 2], Det1[Det1.size() - 1], PulseWidth)) {
							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(d);
							Det1eventValid.push_back(eventIndex);

						}

					}
					else {

						if ((NumberOfCounts1 + (0.5 * NumberOfCounts1_511)) >= (window_511_low)) {

							std::pair<double, double> a = { (Det1_511_arrivalTime + time_fingerprint) , -NumberOfCounts1_511 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det1_1275_arrivalTime + time_fingerprint) , -NumberOfCounts1 }; // time stamp (ns), amplitude (mV) pair
							A.push_back(d);
							Det1eventValid.push_back(eventIndex);
							Det1.push_back(Det1_1275_arrivalTime + time_fingerprint);
							Det1numberOfCounts.push_back(NumberOfCounts1);
						}
					}

					if (Det2.size() > 2) {

						if ((Det2valid[Det2valid.size() - 1] == true) && (Det2valid[Det2valid.size() - 2] == false) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> b = { (Det2[Det2.size() - 2]) , -(Det2numberOfCounts[Det2numberOfCounts.size() - 2]) }; // time stamp (ns), amplitude (mV) pair
							B.push_back(b);
							Det2eventValid.push_back((Det2event[Det2event.size() - 2]));
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(d);
							Det2eventValid.push_back(eventIndex);

						}

						if ((Det2valid[Det2valid.size() - 2] == true) && lookForPileUp(Det2[Det2.size() - 2], Det2[Det2.size() - 1], PulseWidth)) {
							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(d);
							Det2eventValid.push_back(eventIndex);

						}

					}
					else {

						if ((NumberOfCounts2 + (0.5 * NumberOfCounts2_511)) >= (window_511_low)) {

							std::pair<double, double> c = { (Det2_511_arrivalTime + time_fingerprint) , -NumberOfCounts2_511 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(c);
							Det2eventValid.push_back(eventIndex);
							std::pair<double, double> d = { (Det2_1275_arrivalTime + time_fingerprint) , -NumberOfCounts2 }; // time stamp (ns), amplitude (mV) pair
							B.push_back(d);
							Det2eventValid.push_back(eventIndex);

						}

					}

				}
				*/
				if ((Det1.size() >= 2) && (Det1valid[Det1valid.size() - 1] == false) && (Det1valid[Det1valid.size() - 2] == false)) {

					if (lookForClosePileUp(Det1[Det1.size() - 1], Det1[Det1.size() - 2], risingEdge)) {

						double value = Det1numberOfCounts[Det1numberOfCounts.size() - 2] + ((1 / (1 + (Det1[Det1.size() - 1] - Det1[Det1.size() - 2]))) * Det1numberOfCounts[Det1numberOfCounts.size() - 1]);

						if (value <= (window_511_low)) {
							std::pair<double, double> a = { (Det1[Det1.size() - 2]) , -Det1numberOfCounts[Det1numberOfCounts.size() - 2] }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det1eventValid.push_back(eventIndex);
							std::pair<double, double> b = { (Det1[Det1.size() - 1]) , -Det1numberOfCounts[Det1numberOfCounts.size() - 1] }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det1eventValid.push_back(eventIndex);
						}
					}

				}

				if ((Det2.size() >= 2) && (Det2valid[Det2valid.size() - 1] == false) && (Det2valid[Det2valid.size() - 2] == false)) {
					if (lookForClosePileUp(Det2[Det2.size() - 1], Det2[Det2.size() - 2], risingEdge)) {
						double value = Det2numberOfCounts[Det2numberOfCounts.size() - 2] + ((1 / (1 + (Det2[Det2.size() - 1] - Det2[Det2.size() - 2]))) * Det2numberOfCounts[Det2numberOfCounts.size() - 1]);

						if (value <= (window_511_low)) {
							std::pair<double, double> a = { (Det2[Det2.size() - 2]) , -Det2numberOfCounts[Det2numberOfCounts.size() - 2] }; // time stamp (ns), amplitude (mV) pair
							A.push_back(a);
							Det2eventValid.push_back(eventIndex);
							std::pair<double, double> b = { (Det2[Det2.size() - 1]) , -Det2numberOfCounts[Det2numberOfCounts.size() - 1] }; // time stamp (ns), amplitude (mV) pair
							A.push_back(b);
							Det2eventValid.push_back(eventIndex);
						}
					}

				}
				
				time_fingerprint += Source_decay;

				if ((Det1.size() == 0) && (Det2.size() == 0)) {
					A.clear();
					B.clear();

					time_fingerprint = 0.0;
					Det1numberOfCounts.clear();
					Det2numberOfCounts.clear();
					Det1.clear();
					Det2.clear();
					Det1event.clear();
					Det2event.clear();
					sequence.clear();
					Det1valid.clear();
					Det2valid.clear();
					Det1eventValid.clear();
					Det2eventValid.clear();

				}






				Classification_ID = 0;


			}
		}
		std::cout << i << std::endl;
		ifile.clear();
		ifile511.clear();

		ifile.seekg(0, std::ios::beg);
		ifile511.seekg(0, std::ios::beg);


	}
	std::cout << "number of events: " << eventIndex << std::endl;
	myfile << "number of 200ns windows :" << Window200nsCounter << endl;
	myfile << "number of true coincidence :" << TrueCoincidenceCounter << endl;
	myfile << "number of false coincidence :" << FalseCoincidenceCounter << endl;
	myfile << "number of True Pile Up :" << TruePileUpCounter << endl;
	myfile << "number of  False Pile Up :" << FalsePileUpCounter << endl;
	myfile << "Number of solo 511 double detection :" << Double511Detection << endl;
	myfile << "Number of 1275 Backscatter coincidence :" << Backscatter1275PileUpCounter << endl;



	// (D) close file stream ...
	ifile.close();
	ifile511.close();
	myfile.close();

	DPulseStreamManager::sharedInstance()->stopAndSave();

	std::cout << "stream finished!";

	DDELETE_SAFETY(pulseGenerator);

}
