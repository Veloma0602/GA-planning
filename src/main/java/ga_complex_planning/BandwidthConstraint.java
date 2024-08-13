package ga_complex_planning;

/**
 * @ClassName BandwidthConstraint
 * @Description 带宽限制函数
 * @Author fwt
 * @Date 2024/1/24 9:35 下午
 * @Version 1.0
 **/

public class BandwidthConstraint {

    // 定值
    private double lambda;
    private double b_max;

    // 构造函数
    public BandwidthConstraint(Chromosome chromosome) {
        this.lambda = chromosome.getLambda();
        this.b_max = chromosome.getB_max();
    }

    // 计算带宽限制函数
    public double calculateBandwidthPenalty(int[] bandwidths) {
        double bandwidthPenalty = 0.0;
        for (double bandwidth : bandwidths) {
            // 计算每个平台的带宽惩罚项
            bandwidthPenalty += Math.pow((b_max - bandwidth), 2);
        }
        // 计算总的带宽限制函数
        return lambda * bandwidthPenalty;
    }
}