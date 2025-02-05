package pitchfork.models;

import beast.base.core.Input;
import beast.base.evolution.tree.IntervalType;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Binomial;
import pitchfork.models.AbstractBetaSkylineDistribution;

public class BetaSkylineDistributionVersion1 extends AbstractBetaSkylineDistribution {

    public BetaSkylineDistributionVersion1() {
        super();
    }



    @Override
    public double calculateLogP() {
        logP = 0.0;

        double t = 0;
        double N = skylinePopulationsInput.get().getValue(0);
        int seenCoalescentEvents = 0;
        int group = 0;
        int totalCoalecents = 0;  //this variable is needed that N is not updated after the last coalescent events because there is no N value left
        for (int i = 0; i < collapsedTreeIntervals.getIntervalCount(); i++) {
            totalCoalecents += collapsedTreeIntervals.getCoalescentEvents(i);
        }
        int counter = 0;     //this variable counts all the coalescent events and is compared to totalCoalescents

        for (int i = 0; i < collapsedTreeIntervals.getIntervalCount(); i++) {
            // Get interval details
            double dt = collapsedTreeIntervals.getInterval(i);
            int n = collapsedTreeIntervals.getLineageCount(i);

            // Waiting time contribution   correct!
            logP += -betaCoalescentModel.getTotalCoalRate(n)*dt/N;


            // Increment time
            t += dt;

            if (collapsedTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                // Beta-coalescent event contribution

                int k = collapsedTreeIntervals.getCoalescentEvents(i)+1;
                logP += betaCoalescentModel.getLogLambda(n, k) + Binomial.logChoose(n, k) - Math.log(N);

                // while loop to update N until the seen coalescent events are smaller or equal than the actual group size
                seenCoalescentEvents += collapsedTreeIntervals.getCoalescentEvents(i);
                counter += collapsedTreeIntervals.getCoalescentEvents(i);

                if (counter < totalCoalecents) {           //because we don't want to update N after the last coalescent event
                    while (seenCoalescentEvents >= groupSizesInput.get().getValue(group)) {
                        seenCoalescentEvents -= groupSizesInput.get().getValue(group);
                        group += 1;
                    }
                    N = skylinePopulationsInput.get().getValue(group);
                }
            }
        }
        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    double[] getPopSizes(int gridSize) {
        double[] popSizes = new double[gridSize];
        double T = tree.getRoot().getHeight();

        int group = 0;
        int seenCoalescentEvents = 0;
        double N = skylinePopulationsInput.get().getValue(0);
        int i = 0;
        double t_total = 0;
        int counter = 0;
        int totalCoalecents = 0;
        for (int j = 0; j < collapsedTreeIntervals.getIntervalCount(); j++) {
            totalCoalecents += collapsedTreeIntervals.getCoalescentEvents(j);
        }

        for (int gridIdx = 0; gridIdx < gridSize; gridIdx++) {
            double t = gridIdx * T / (gridSize - 1);
            //int coalescent = 0;

            // TODO: update N as needed
            while((t_total+collapsedTreeIntervals.getInterval(i)) < t){

                double dt = collapsedTreeIntervals.getInterval(i);

                // Increment time
                t_total += dt;

                if (collapsedTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    seenCoalescentEvents += collapsedTreeIntervals.getCoalescentEvents(i);
                    counter += collapsedTreeIntervals.getCoalescentEvents(i);
                    if (counter < totalCoalecents) {
                        while (seenCoalescentEvents > groupSizesInput.get().getValue(group)) {    // > not >= because above N is for next round but not here
                            seenCoalescentEvents -= groupSizesInput.get().getValue(group);
                            group += 1;
                        }
                        N = skylinePopulationsInput.get().getValue(group);
                    }
                }

                if (i == (collapsedTreeIntervals.getIntervalCount()-1)) {  //this if statement ends while loop when we reach the root of the tree
                    break;
                }
                i++;
            }

            popSizes[gridIdx] = N;
        }

        return popSizes;
    }
}
