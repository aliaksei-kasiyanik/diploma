package com.akasiyanik.bsu.coursework.methods;

import com.akasiyanik.bsu.coursework.equations.SteadyingEquation;
import com.akasiyanik.bsu.coursework.methods.rungekutta.AuxRungeKuttaMethod;

/**
 * @author: akasiyanik
 */
public class AdaptStepsizeSteadyingProcess extends SteadyingProcess {

    public AdaptStepsizeSteadyingProcess(AuxRungeKuttaMethod auxMethod,
                                         double EPS,
                                         double w,
                                         SteadyingEquation steadyingEquation) {
        super(auxMethod, EPS, w, steadyingEquation);
    }
}
