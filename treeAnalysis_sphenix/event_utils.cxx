/**
  * Event level utilities
  *
  * - Contains things like DIS kinematics
  *
  * \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
  */

struct DISKinematics {
    DISKinematics() = default;

    double x;
    double y;
    double Q2;
};

/**
  * Determine DIS kinematics variables using the Jacquet-Blondel method
  *
  * Implemented using https://wiki.bnl.gov/eic/index.php/DIS_Kinematics
  */
DISKinematics JBKinematics(unsigned short primaryTrackSource, double incomingElectronEnergy, double incomingHadronEnergy, double hadronMassHypothesis = 0.139)
{
    // NOTE: We need electron ID here. Since we're don't yet have that available as of 16 June 2021, we use the true PID info.
    // NOTE: 1e6 is a sentinel value.
    unsigned int electronID = 1e6;
    double electronPt = -1e6;
    for (unsigned int i = 0; i < (unsigned int)m_nTracks; ++i) {

        if (_mcpart_PDG[static_cast<int>(m_truthtrackID[i])] == 11) {
            // If it's an electron, take the highest pt one.
            // TODO: Update this when the electron PID is addressed.
            double pt = std::sqrt(std::pow(m_tr_px[i], 2) + std::pow(m_tr_py[i], 2));
            if (pt > electronPt) {
                electronID = i;
                electronPt = pt;
            }
        }
    }
    if (electronID == 1e6) {
        throw std::runtime_error("Unable to find scattered electron!");
    }

    // Loop over the final state hadrons to determine the kinematics.
    double px_h = 0;
    double py_h = 0;
    double e_minus_pz = 0;
    for(unsigned int i = 0; i < (unsigned int)m_nTracks; ++i) {
        // Only consider tracks from the primary track source.

        // Skip the scattered electron.
        if (i == electronID) {
            continue;
        }
        px_h += m_tr_px[i];
        py_h += m_tr_py[i];
        double eTrack = std::sqrt(std::pow(m_tr_px[i],2) + std::pow(m_tr_py[i],2) + std::pow(m_tr_pz[i],2) + std::pow(hadronMassHypothesis, 2));
        e_minus_pz += (eTrack - m_tr_pz[i]);
    }

    // TODO: Add calorimeter towers.

    double sqrts = 4 * incomingElectronEnergy * incomingHadronEnergy;

    double y = e_minus_pz / (2 * incomingElectronEnergy);
    double Q2 = (std::pow(px_h, 2) + std::pow(py_h, 2)) / (1 - y);
    double x = Q2 / (sqrts * y);

    return DISKinematics{x, y, Q2};
}
