#pragma once
#ifndef DETEVENTINFO_HH
#define DETEVENTINFO_HH

using namespace std;

double positronAnhillationLifetime(int NumberOfComponents, double lifetime1, double lifetime2, double lifetime3, double lifetime4, double lifetime5, double intensity1,
    double intensity2, double intensity3, double intensity4, double intensity5, std::exponential_distribution<double>& distribution1,
    std::exponential_distribution<double>& distribution2, std::exponential_distribution<double>& distribution3,
    std::exponential_distribution<double>& distribution4, std::exponential_distribution<double>& distribution5,
    std::mt19937& generator)

{   double positronLifetime;
    
if (NumberOfComponents == 1) {
    double overall_intensity = intensity1;
    double tolerance = 1e-9; // set a small tolerance value
    if (std::abs(overall_intensity - 1.0) > tolerance) {
        throw std::invalid_argument("Error: The intensities must sum up to 1!");
    }


    double r = ((double)rand() / RAND_MAX);
    if (r <= intensity1 / overall_intensity)
    {
        distribution1.param(std::exponential_distribution<double>::param_type(1 / lifetime1));

        positronLifetime = (distribution1(generator));
    }
    }else if (NumberOfComponents == 2) {
        double overall_intensity = intensity1 + intensity2 + intensity3;
        double tolerance = 1e-9; // set a small tolerance value
        if (std::abs(overall_intensity - 1.0) > tolerance) {
            throw std::invalid_argument("Error: The intensities must sum up to 1!");
        }


        double r = ((double)rand() / RAND_MAX);
        if (r <= intensity1 / overall_intensity)
        {
            distribution1.param(std::exponential_distribution<double>::param_type(1 / lifetime1));

            positronLifetime = (distribution1(generator));
        }
        else if (r <= (intensity1 + intensity2) / overall_intensity)
        {
            distribution2.param(std::exponential_distribution<double>::param_type(1 / lifetime2));

            positronLifetime = (distribution2(generator));
        }
        else {
            distribution3.param(std::exponential_distribution<double>::param_type(1 / lifetime3));
            positronLifetime = (distribution3(generator));
        }
    }else if(NumberOfComponents == 3) {
        double overall_intensity = intensity1 + intensity2 + intensity3;
        double tolerance = 1e-9; // set a small tolerance value
        if (std::abs(overall_intensity - 1.0) > tolerance) {
            throw std::invalid_argument("Error: The intensities must sum up to 1!");
        }

        
        double r = ((double)rand() / RAND_MAX);
        if (r <= intensity1 / overall_intensity)
        {
            distribution1.param(std::exponential_distribution<double>::param_type(1 / lifetime1));

            positronLifetime = (distribution1(generator));
        }
        else if (r <= (intensity1 + intensity2) / overall_intensity)
        {
            distribution2.param(std::exponential_distribution<double>::param_type(1 / lifetime2));

            positronLifetime = (distribution2(generator));
        }
        else {
            distribution3.param(std::exponential_distribution<double>::param_type(1 / lifetime3));
            positronLifetime = (distribution3(generator));
        }
    }else if (NumberOfComponents == 4) {
        double overall_intensity = intensity1 + intensity2 + intensity3 + intensity4;
        double tolerance = 1e-9; // set a small tolerance value
        if (std::abs(overall_intensity - 1.0) > tolerance) {
            throw std::invalid_argument("Error: The intensities must sum up to 1!");
        }


        double r = ((double)rand() / RAND_MAX);
        if (r <= intensity1 / overall_intensity)
        {
            distribution1.param(std::exponential_distribution<double>::param_type(1 / lifetime1));

            positronLifetime = (distribution1(generator));
        }
        else if (r <= (intensity1 + intensity2) / overall_intensity)
        {
            distribution2.param(std::exponential_distribution<double>::param_type(1 / lifetime2));

            positronLifetime = (distribution2(generator));
        }
        else if (r <= (intensity1 + intensity2 + intensity3) / overall_intensity)
        {
            distribution3.param(std::exponential_distribution<double>::param_type(1 / lifetime3));

            positronLifetime = (distribution3(generator));
        }
        else {
            distribution4.param(std::exponential_distribution<double>::param_type(1 / lifetime4));
            positronLifetime = (distribution4(generator));
        }
    }else if (NumberOfComponents == 4) {
        double overall_intensity = intensity1 + intensity2 + intensity3 + intensity4 + intensity5;
        double tolerance = 1e-9; // set a small tolerance value
        if (std::abs(overall_intensity - 1.0) > tolerance) {
            throw std::invalid_argument("Error: The intensities must sum up to 1!");
        }


        double r = ((double)rand() / RAND_MAX);
        if (r <= intensity1 / overall_intensity)
        {
            distribution1.param(std::exponential_distribution<double>::param_type(1 / lifetime1));

            positronLifetime = (distribution1(generator));
        }
        else if (r <= (intensity1 + intensity2) / overall_intensity)
        {
            distribution2.param(std::exponential_distribution<double>::param_type(1 / lifetime2));

            positronLifetime = (distribution2(generator));
        }
        else if (r <= (intensity1 + intensity2 + intensity3) / overall_intensity)
        {
            distribution3.param(std::exponential_distribution<double>::param_type(1 / lifetime3));

            positronLifetime = (distribution3(generator));
        }
        else if (r <= (intensity1 + intensity2 + intensity3 + intensity4) / overall_intensity)
        {
            distribution4.param(std::exponential_distribution<double>::param_type(1 / lifetime4));

            positronLifetime = (distribution4(generator));
        }
        else {
            distribution5.param(std::exponential_distribution<double>::param_type(1 / lifetime5));
            positronLifetime = (distribution5(generator));
        }
        }
    return positronLifetime;
}

double generatePMTuncertaincy(double timeTransitSpread, double QE, int numberOfCounts) {
    if (numberOfCounts <= 3) {
        return 0.0;
    }
    double value = timeTransitSpread * (1/((sqrt(numberOfCounts * QE))));
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::normal_distribution<double> distribution((0.0, value));

    if (distribution.param() != std::normal_distribution<double>().param()) {
        distribution.param(std::normal_distribution<double>::param_type(0.0, value));
    }
    
    return distribution(gen);
}

double radioactivity(double decay_constant, std::exponential_distribution<double>& distribution, std::mt19937& generator) {

    distribution.param(std::exponential_distribution<double>::param_type(1 / decay_constant));

    double DecayTime = (distribution(generator)) * 1e9;

    return DecayTime;
}

class DetectorInfo { // A || B
    friend class EventInfo;
public:
    DetectorInfo() {
        clear();
    }

    DetectorInfo(int id) {
        clear();

        m_id = id;
    }

    ~DetectorInfo() {}

    void clear() {
        m_id = -1;
        m_depEnerg = 0;
        //m_ArrCF = 0.;
        m_MeanArr = 0.;
        //m_MedianArr = 0.;
       // m_MaxPulse = 0.;
        m_globalTime = 0;
        m_nofHits = false;
        m_interactionX = 0.;
        m_interactionY = 0.;
        m_interactionZ = 0.;
        m_numberOfCounts = 0;
        m_MeanTrackLength = 0;

    }


    void AddEnergyDep(double fEdep, double globalTime, bool bHit, int id, double MeanArr, double InterActionX, double InterActionY,
        double InterActionZ, int numberOfCounts, double MeanTrackLength, int NumberOfReflection, int GammaAgain/*, double MedianArr, double MaxPulse*/) {
        m_depEnerg = fEdep;
        //m_ArrCF = ArrCF;
        m_globalTime = globalTime;
        m_nofHits = bHit;
        m_id = id;
        m_MeanArr = MeanArr;
        m_interactionX = InterActionX;
        m_interactionY = InterActionY;
        m_interactionZ = InterActionZ;
        m_numberOfCounts = numberOfCounts;
        m_MeanTrackLength = MeanTrackLength;
        m_NumberOfReflection = NumberOfReflection;
        m_GammaAgain = GammaAgain;

        // m_MedianArr = MedianArr;
        // m_MaxPulse = MaxPulse;
    }


    int NumberOfReflection() const {
        return m_NumberOfReflection;
    }
    int GammaAgain() const {
        return m_GammaAgain;
    }

    int id() const {
        return m_id;
    }

    double energy() const {
        return m_depEnerg;
    }

    /* double ArrCF() const {
         return m_ArrCF;
     }*/

    double MeanArr() const {
        return m_MeanArr;
    }
    double InterActionX() const {
        return m_interactionX;
    }
    double InterActionY() const {
        return m_interactionY;
    }
    double InterActionZ() const {
        return m_interactionZ;
    }
    int numberOfCounts() const {
        return m_numberOfCounts;
    }

    /*
         double MedianArr() const {
            return m_MedianArr;
        }

        double MaxPulse() const {
            return m_MaxPulse;
        }
    */

    double MeanTrackLength() const {
        return m_MeanTrackLength;
    }

    bool isValid() const {
        return m_nofHits;
    }

private:
    int m_id; // event-id
    double m_depEnerg; // dep. energy at detector
    //double m_ArrCF; // arrival time at detector
    double m_MeanArr; //mean arrvialValue
    double m_globalTime; //globalTime counter
    double m_MedianArr; //Median value of the pulse
    double m_MaxPulse; //Pulse Peak x value
    bool m_nofHits; // incident gammay ray?
    double m_interactionX;
    double m_interactionY;
    double m_interactionZ;
    int m_numberOfCounts;

    double m_MeanTrackLength;

    int m_NumberOfReflection;
    int m_GammaAgain;


};

class EventInfo {
public:
    EventInfo() { clear(); }
    ~EventInfo() {}

    void clear() {
        m_lifetime = 0.;
        m_startTime = 0.;
        m_StartWinkelX = 0.;
        m_StartWinkelY = 0.;
        m_StartWinkelZ = 0.;
        m_StoppWinkel1X = 0.;
        m_StoppWinkel1Y = 0.;
        m_StoppWinkel1Z = 0.;
        m_StoppWinkel2X = 0.;
        m_StoppWinkel2Y = 0.;
        m_StoppWinkel2Z = 0.;
        m_detector1.clear();
        m_detector2.clear();
    }

    void setInfo(double lifetime, double startTime, float StartWinkelX, float StartWinkelY, float StartWinkelZ, float StoppWinkel1X,
        float StoppWinkel1Y, float StoppWinkel1Z, float StoppWinkel2X, float StoppWinkel2Y, float StoppWinkel2Z) {
        m_lifetime = lifetime;
        m_startTime = startTime;
        m_StartWinkelX = StartWinkelX;
        m_StartWinkelY = StartWinkelY;
        m_StartWinkelZ = StartWinkelZ;
        m_StoppWinkel1X = StoppWinkel1X;
        m_StoppWinkel1Y = StoppWinkel1Y;
        m_StoppWinkel1Z = StoppWinkel1Z;
        m_StoppWinkel2X = StoppWinkel2X;
        m_StoppWinkel2Y = StoppWinkel2Y;
        m_StoppWinkel2Z = StoppWinkel2Z;
    }

    void attach(const DetectorInfo& detector1, const DetectorInfo& detector2) {
        m_detector1 = detector1;
        m_detector2 = detector2;
    }

    DetectorInfo info1() const {
        return m_detector1;
    }

    DetectorInfo info2() const {
        return m_detector2;
    }


    double lifetime() const {
        return m_lifetime;
    }

    double startTime() const {
        return m_startTime;
    }

    float StartWinkelX() const {
        return m_StartWinkelX;
    }

    float StartWinkelY() const {
        return m_StartWinkelY;
    }

    float StartWinkelZ() const {
        return m_StartWinkelZ;
    }

    float StoppWinkel1X() const {
        return m_StoppWinkel1X;
    }
    float StoppWinkel1Y() const {
        return m_StoppWinkel1Y;
    }
    float StoppWinkel1Z() const {
        return m_StoppWinkel1Z;
    }

    float StoppWinkel2X() const {
        return m_StoppWinkel2X;
    }
    float StoppWinkel2Y() const {
        return m_StoppWinkel2Y;
    }
    float StoppWinkel2Z() const {
        return m_StoppWinkel2Z;
    }

private:
    DetectorInfo m_detector1;
    DetectorInfo m_detector2;

    double m_lifetime;
    double m_startTime;
    float m_StartWinkelX;
    float m_StartWinkelY;
    float m_StartWinkelZ;
    float m_StoppWinkel1X;
    float m_StoppWinkel1Y;
    float m_StoppWinkel1Z;
    float m_StoppWinkel2X;
    float m_StoppWinkel2Y;
    float m_StoppWinkel2Z;
};

#endif
