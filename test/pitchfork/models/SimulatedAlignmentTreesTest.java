package pitchfork.models;

import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import pitchfork.models.BetaCoalescentDistribution;
import feast.simulation.SimulatedAlignment;
import org.junit.Test;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.substitutionmodel.SubstitutionModel;


public class SimulatedAlignmentTreesTest {

    @Test
    public void testAlignment() {
        TreeParser tree = new TreeParser("((A:1.0,B:1.0,C:1.0):0.5,D:0.7);",
                false, false, true, 0);
        SimulatedAlignment saModel = new SimulatedAlignment();
        JukesCantor jc = new JukesCantor();
        saModel.initByName("tree", tree,
                "siteModel", jc,
                "sequenceLength", new IntegerParameter("15"),
                "outputFileName", "SimulatedAlignment.trees");
    }

}
