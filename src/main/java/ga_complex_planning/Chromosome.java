package ga_complex_planning;


import ga_complex_planning.planning_info.Info;
import ga_complex_planning.pojo.Point;
import util.GAGraphUtil;
import util.PropertyUtil;
import java.util.*;


/**
 * @ClassName Chromosome
 * @Description TODO
 * @Author fwt
 * @Date 2024/1/24 9:35 下午
 * @Version 1.0
 **/
public class Chromosome {
    private int[] gene;
    private int geneSize;
    private double score;
    private double[] frequencies; // 移动平台的频率数组
    private double[] bandwidths;  // 移动平台的带宽数组
    // 将参数定义为类的成员变量
    private int numOfPlatforms;

    private double f_min;  //频率上下限
    private double f_max;
    private double f_disable_min;     //禁用频率上下限
    private double f_disable_max;
    private double b_range_min;       //带宽上下限
    private double b_range_max;
    private double b_max;
    private double fixedTotalPowerValue;    //总功率

    

    // 带宽限制函数中参数值
    private double lambda = 2.0;
    

    


    private static Random r = new Random();

    public Chromosome() {}


    // 构造函数接收这些参数
    public Chromosome(int numOfPlatforms, double f_min, double f_max, 
                       double f_disable_min, double f_disable_max,
                       double b_range_min, double b_range_max, double b_max
                       ) {
        this.numOfPlatforms = numOfPlatforms;
        this.f_min = f_min;
        this.f_max = f_max;
        this.f_disable_min = f_disable_min;
        this.f_disable_max = f_disable_max;
        this.b_range_min = b_range_min;
        this.b_range_max = b_range_max;
        this.b_max = b_max;
        frequencies = new double[numOfPlatforms];
        bandwidths = new double[numOfPlatforms];
        init();

    }


    /**
     * 初始化个体基因，需要保证
     * 1、基因中起始点0个数 为 carNum +1
     * 2、每个点只能出现一次，且所有点都必须出现
     * 3、两个起始点不可以相邻
     * 4、基因收尾必须是起始点位置 0
     *
     * 假设 carNum = 2  pointNum =4
     * 正确的基因:
     * 0 2 3 0 1 4 0
     *
     * 错误的基因：
     * 0 3 3 0 1 4 0  (3受灾点重复)
     * 0 0 3 2 1 4 0  (有两个起始点重合->有车子没出发)
     */
        
        public void init() {
            Random random = new Random();
            double totalBandwidth = 0.0;
    
            for (int i = 0; i < frequencies.length; i++) {
                frequencies[i] = f_min + random.nextDouble() * (f_max - f_min);
    
                // 确保频率不在禁用频段内
                while (frequencies[i] >= f_disable_min && frequencies[i] <= f_disable_max) {
                    frequencies[i] = f_min + random.nextDouble() * (f_max - f_min);
                }
    
                bandwidths[i] = b_range_min + random.nextDouble() * (b_range_max - b_range_min);
                totalBandwidth += bandwidths[i];

            }
    
            // 确保带宽总和不超过可分配总带宽
//            double totalBandwidth = Arrays.stream(bandwidths).sum();

            double scale = Math.min(1.0, b_max / totalBandwidth - 0.05); // 计算调整比例，确保不超过上限
            double DtotalBandwidth = 0.0;
            for (int i = 0; i < bandwidths.length; i++) {
                bandwidths[i] *= scale;
                DtotalBandwidth += bandwidths[i];
//                System.out.println(bandwidths[i]);
            }
//            System.out.println("总带宽为"+DtotalBandwidth);
        }
    


    /**
     * 判断当前染色体是否符合要求
     * 1、基因中起始点0个数 为 carNum +1
     * 2、每个点只能出现一次，且所有点都必须出现
     * 3、两个起始点不可以相邻
     * 4、基因收尾必须是起始点位置 0
     * 这个函数主要是用在交叉后对子代进行判断的
     * 这里还要注意一点，如果生成的子代不满足条件
     * 我们不应该重新轮盘赌选择父代，而是重新选择交叉片段的位移
     * @param chromosome
     * @return
     */
  

    /**
     * 交叉生成两个新的子代
     * 这里确保两个子代都是符合染色体规范的
     * @param p1
     * @param p2
     * @return
     */

     public static List<Chromosome> genetic(Chromosome parent1, Chromosome parent2) {
        if (parent1 == null || parent2 == null || parent1.getNumOfPlatforms() != parent2.getNumOfPlatforms()) {
            return null;
        }

        int numOfPlatforms = parent1.getNumOfPlatforms();
        double f_min = parent1.getF_min();
        double f_max = parent1.getF_max();
        double f_disable_min = parent1.getF_disable_min();
        double f_disable_max = parent1.getF_disable_max();
        double b_range_min = parent1.getB_range_min();
        double b_range_max = parent1.getB_range_max();
        double b_max = parent1.getB_max();

        int crossoverPoint = (int) (Math.random() * parent1.getNumOfPlatforms()); // 交叉点

        Chromosome child1 = new Chromosome(numOfPlatforms, f_min, f_max, f_disable_min, f_disable_max,
                                            b_range_min, b_range_max, b_max);
        Chromosome child2 = new Chromosome(numOfPlatforms, f_min, f_max, f_disable_min, f_disable_max,
                                            b_range_min, b_range_max, b_max);

        double[] parent1Frequencies = parent1.getFrequencies();
        double[] parent2Frequencies = parent2.getFrequencies();
        double[] parent1Bandwidths = parent1.getBandwidths();
        double[] parent2Bandwidths = parent2.getBandwidths();

        // 交叉频率信息
        for (int i = 0; i < crossoverPoint; i++) {
            child1.setFrequency(i, parent1Frequencies[i]);
            child2.setFrequency(i, parent2Frequencies[i]);
        }
        for (int i = crossoverPoint; i < numOfPlatforms; i++) {
            child1.setFrequency(i, parent2Frequencies[i]);
            child2.setFrequency(i, parent1Frequencies[i]);
        }

        // 交叉带宽信息
        for (int i = 0; i < crossoverPoint; i++) {
            child1.setBandwidth(i, parent1Bandwidths[i]);
            child2.setBandwidth(i, parent2Bandwidths[i]);
        }
        for (int i = crossoverPoint; i < numOfPlatforms; i++) {
            child1.setBandwidth(i, parent2Bandwidths[i]);
            child2.setBandwidth(i, parent1Bandwidths[i]);
        }

        List<Chromosome> children = new ArrayList<>();
        children.add(child1);
        children.add(child2);
        return children;
    }


    /**
     * 基因变异
     * 为保证单一变量使得适应度因解决方案而改变
     * 基因内只有带宽发生突变，频率不发生突变
     *
     * @param maxMutationNum
     */
    // TODO:
    public void mutation(int maxMutationNum,Chromosome chromosome) {
        // 如果要变异，变异对的个数起码也得是1对
        int mutationPairNum = Math.max((int) Math.random() * maxMutationNum / 2,1);
        // 变异策略有两种
        // 1、找出所有合理变异对，然后选择
        // 2、随机选择变异对，然后判断是否合理
        // 这里选择第二种，因为其命中概率不低，不会造成太高的迭代

        for (int i = 0; i < mutationPairNum; i++) {
            int left = r.nextInt(numOfPlatforms - 2) + 1;
            int right = r.nextInt(numOfPlatforms - 2) + 1;
            // 防止出现过高迭代，导致阻塞
            int maxIterNum = 400;
            while (!isSwapPairLegal(chromosome,left,right) && maxIterNum>=0) {
                left = r.nextInt(numOfPlatforms - 2) + 1;
                right = r.nextInt(numOfPlatforms - 2) + 1;
                maxIterNum--;
                if (maxIterNum<0) {
                    left=0;
                    right=0;
                }
            }
            // 交换目标基因
            double Btmp = chromosome.bandwidths[left];
            chromosome.bandwidths[left] = chromosome.bandwidths[right];
            chromosome.bandwidths[right] = Btmp;
        }
    }

    /**
     * 判断待变异前的数对是否符合要求
     * @param left
     * @param right
     * @return
     */
    private boolean isSwapPairLegal(Chromosome chromosome,int left,int right) {
        // 都是非 startPoint，可以变异
        if (chromosome.bandwidths[left]!=0 && chromosome.bandwidths[right]!=0) return true;
        // 无效变异
        if (chromosome.bandwidths[left]==0 && chromosome.bandwidths[right]==0) return false;
        // // 判断交换后是否还合法
        // Chromosome dummy = this.clone();
        // int[] gene = dummy.getGene();
        // int tmp = gene[left];
        // gene[left]=gene[right];
        // gene[right]=tmp;
        // return isGoodChromosome(dummy);
        return true;
    }

    /**
     * 对染色体进行深拷贝
     * @return
     */
    // public Chromosome clone() {
    //     Chromosome child = new Chromosome();
    //     child.setCarNum(carNum);
    //     child.setPointNum(pointNum);
    //     child.setGeneSize(geneSize);
    //     int[] pGene = gene;
    //     int[] cGene = new int[pGene.length];
    //     for (int i = 0; i < pGene.length; i++) {
    //         cGene[i]=pGene[i];
    //     }
    //     child.setGene(cGene);
    //     return child;
    // }

    public int[] getGene() {
        return gene;
    }

    public void setGene(int[] gene) {
        this.gene = gene;
    }

    public int getGeneSize() {
        return geneSize;
    }

    public void setGeneSize(int geneSize) {
        this.geneSize = geneSize;
    }

    public double getScore() {
        return score;
    }

    public void setScore(double score) {
        this.score = score;
    }

    public int getNumOfPlatforms() {
        return numOfPlatforms;
    }
    
    public void setNumOfPlatforms(int numOfPlatforms) {
        this.numOfPlatforms = numOfPlatforms;
    }
    
    public double getF_min() {
        return f_min;
    }
    
    public void setF_min(double f_min) {
        this.f_min = f_min;
    }
    
    public double getF_max() {
        return f_max;
    }
    
    public void setF_max(double f_max) {
        this.f_max = f_max;
    }
    
    public double getF_disable_min() {
        return f_disable_min;
    }
    
    public void setF_disable_min(double f_disable_min) {
        this.f_disable_min = f_disable_min;
    }
    
    public double getF_disable_max() {
        return f_disable_max;
    }
    
    public void setF_disable_max(double f_disable_max) {
        this.f_disable_max = f_disable_max;
    }
    
    public double getB_range_min() {
        return b_range_min;
    }
    
    public void setB_range_min(double b_range_min) {
        this.b_range_min = b_range_min;
    }
    
    public double getB_range_max() {
        return b_range_max;
    }
    
    public void setB_range_max(double b_range_max) {
        this.b_range_max = b_range_max;
    }
    
    public double getB_max() {
        return b_max;
    }
    
    public void setB_max(double b_max) {
        this.b_max = b_max;
    }
    
    public double getFixedTotalPowerValue() {
        return fixedTotalPowerValue;
    }
    
    public void setFixedTotalPowerValue(double fixedTotalPowerValue) {
        this.fixedTotalPowerValue = fixedTotalPowerValue;
    }
    
    public double[] getFrequencies() {
        return frequencies;
    }
    
    public void setFrequencies(double[] frequencies) {
        this.frequencies = frequencies;
    }
    
    public double[] getBandwidths() {
        return bandwidths;
    }
    
    public void setBandwidths(double[] bandwidths) {
        this.bandwidths = bandwidths;
    }


    public double getLambda() {
        return lambda;
    }

    public void setLambda(double lambda) {
        this.lambda = lambda;
    }

    public void setFrequency(int index, double frequency) {
        frequencies[index] = frequency;
    }
    public void setBandwidth(int index, double bandwidth) {
        bandwidths[index] = bandwidth;
    }



    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        return builder.append("基因为:").append(Arrays.toString(gene)).append("\n")
                .append("适应度为:").append(score)
                .toString();
    }
}
