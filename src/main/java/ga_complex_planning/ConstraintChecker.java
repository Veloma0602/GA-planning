package ga_complex_planning;

/**
 * @ClassName ConstraintChecker
 * @Description 检查基因组是否超出约束
 * @Author fwt
 * @Date 2024/1/25 4:35 下午
 * @Version 1.0
 **/
public class ConstraintChecker {
    private static int i = 0;

    public static boolean isGoodChromosome(Chromosome chromosome) {
        // 频段约束
        if (!checkFrequencyRangeConstraint(chromosome.getFrequencies(), chromosome.getF_min(), chromosome.getF_max())) {
            System.out.println("频率出界");
            return false;
        }

        // 频点范围约束
        if (!checkFrequencyPointRangeConstraint(chromosome.getFrequencies(), chromosome.getF_disable_min(), chromosome.getF_disable_max())) {
            System.out.println("频率禁忌");
            return false;

        }

        // 带宽上限约束
        if (!checkBandwidthLimitConstraint(chromosome.getBandwidths(), chromosome.getB_max())) {
//            System.out.println("带宽出界");
            return false;
        }

        return true;
    }

    private static boolean checkFrequencyRangeConstraint(double[] frequencies, double fMin, double fMax) {
        for (double frequency : frequencies) {
            if (frequency < fMin || frequency > fMax) {

                return false;
            }
        }
        return true;
    }

    private static boolean checkFrequencyPointRangeConstraint(double[] frequencies, double fRangeMin, double fRangeMax) {
        for (double frequency : frequencies) {
            if (frequency >= fRangeMin && frequency <= fRangeMax) {

                return false; // 如果有任何一个频率在范围内，返回 false，表示不符合约束
            }
        }
        return true; // 所有频率都不在范围内，返回 true，表示满足约束
    }

    private static boolean checkBandwidthLimitConstraint(double[] bandwidths, double bMax) {
        i++;
        double totalBandwidth = 0.0;
        for (double bandwidth : bandwidths) {
//            System.out.println(bandwidth);
            totalBandwidth += bandwidth;
        }
//        System.out.println("第"+i+"个基因总带宽为"+totalBandwidth);
        return totalBandwidth <= bMax;
    }

}
