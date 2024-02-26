package ga_complex_planning;

import ga_complex_planning.planning_info.Info;
import ga_complex_planning.pojo.Point;
import util.GAGraphUtil;
import util.PropertyUtil;
import java.util.*;



/**
 * @ClassName GeneticAlgorithm
 * @Description 遗传算法
 * @Author fwt
 * @Date 2024/1/25 4:35 下午
 * @Version 1.0
 **/
public class GeneticAlgorithm   {

    private Properties gaComplexPro = PropertyUtil.getProperty("ga_complex");
    private Properties planningInfoPro = PropertyUtil.getProperty("planning_info");

    private Map<Double,Double> bestScoreDataMap = new HashMap<>();
    private Map<Double,Double> worstScoreDataMap = new HashMap<>();
    private Map<Double,Double> totalScoreDataMap = new HashMap<>();

    public static int underZeroCount = 0;

    private static final Random r = new Random();


    // 受灾地信息
    private Map<String, Point> info = Info.getInfo();

    double wInterference = Double.valueOf(gaComplexPro.getProperty("wInterference"));
    double wBandwidth = Double.valueOf(gaComplexPro.getProperty("wBandwidth")); //权重还没设置
    
    double ORIGIN_SCORE = Double.valueOf(gaComplexPro.getProperty("ORIGIN_SCORE"));
    // 种群
    private List<Chromosome> pop = new ArrayList<>();
    // 种群大小
    private int POP_SIZE = Integer.valueOf(gaComplexPro.getProperty("POP_SIZE"));
    // 迭代次数
    private int ITER_NUM = Integer.valueOf(gaComplexPro.getProperty("ITER_NUM"));
    // 迭代次数计数
    private int iterCount = 0;

    private int numOfPlatforms = Integer.valueOf(planningInfoPro.getProperty("numOfPlatforms"));

    private double f_min = Double.valueOf(planningInfoPro.getProperty("f_min"));  //频率上下限
    private double f_max = Double.valueOf(planningInfoPro.getProperty("f_max"));  //频率上下限
    private double f_disable_min = Double.valueOf(planningInfoPro.getProperty("f_disable_min"));     //禁用频率上下限
    private double f_disable_max = Double.valueOf(planningInfoPro.getProperty("f_disable_max"));
    private double b_range_min = Double.valueOf(planningInfoPro.getProperty("b_range_min"));       //带宽上下限
    private double b_range_max = Double.valueOf(planningInfoPro.getProperty("b_range_max"));
    private double b_max = Double.valueOf(planningInfoPro.getProperty("b_max"));
    private double fixedTotalPowerValue = Double.valueOf(planningInfoPro.getProperty("fixedTotalPowerValue"));    //总功率


    // 当前种群最佳适应度
    private double bestScore = Double.MIN_VALUE;
    // 当前种群最坏适应度
    private double worstScore = Double.MAX_VALUE;
    // 当前种群总适应度
    private double totalScore = 0;
    // 当前种群平均适应度
    private double averageScore = 0;
    // 当前种群最佳基因
    private Chromosome bestGene = null;
    // 当前种群最坏基因
    private Chromosome worstGene = null;

    // 变异概率
    private double MUTATION_RATE = Double.valueOf(gaComplexPro.getProperty("MUTATION_RATE"));
    // 最大变异长度
    private int MAX_MUTATION_NUM = 2;


    /**
     * 执行遗传算法函数
     * @return 执行结束后最佳的种群基因
     */
    public Chromosome conductGA() {
        long startTime = System.currentTimeMillis();
        init();
        for (int i = 0; i < ITER_NUM; i++) {
            // 1、计算种群适应度
            calculatePopScore();
            print(i+1);
            // 2、交叉生成新的种群
            evolve();
             // 3、种群变异
             mutation();
        }
        long endTime = System.currentTimeMillis();
        System.out.println("==============================");
        System.out.println("出现负值结果个数为："+GeneticAlgorithm.underZeroCount);
        System.out.println("平均每次迭代出现负值的个数为:"+GeneticAlgorithm.underZeroCount/500);
        System.out.println("运算耗时为："+(endTime-startTime)+"ms");
        Map[] bestWorstDataSet = new Map[2];
        Map[] totalDataSet = new Map[1];
        bestWorstDataSet[0]=bestScoreDataMap;
        bestWorstDataSet[1]=worstScoreDataMap;
        totalDataSet[0]=totalScoreDataMap;
        GAGraphUtil.drawBestWorstScoreGraph(bestWorstDataSet);
        GAGraphUtil.drawTotalScoreGraph(totalDataSet);
        GAGraphUtil.blockUtil();
        return bestGene;
    }

    /**
     * 临时打印相关信息
     */
    private void print(int generation) {
        System.out.println("----------------------------------------");
        System.out.println("当前代数："+generation);
        System.out.println("当前种群最佳适应度："+bestScore);
        System.out.println("当前种群最坏适应度："+worstScore);
        System.out.println("当前种群总适应度："+totalScore);
        System.out.println("当前种群平均适应度："+averageScore);
        System.out.println("当前种群最佳基因的频率分布："+Arrays.toString(bestGene.getFrequencies()));
        System.out.println("当前种群最佳基因的带宽分布："+Arrays.toString(bestGene.getBandwidths()));
        System.out.println("当前种群最差基因的频率分布："+Arrays.toString(worstGene.getFrequencies()));
        System.out.println("当前种群最差基因的带宽分布："+Arrays.toString(worstGene.getBandwidths()));

    }

    /**
     * 初始化种群
     */
    private void init() {
        for (int i = 0; i < POP_SIZE ;i++) {
            Chromosome chromosome = new Chromosome(numOfPlatforms,f_min,f_max, 
            f_disable_min,f_disable_max,
            b_range_min,b_range_max,b_max);
            pop.add(chromosome);
        }
    }

    /**
     * 计算种群适应度
     */
    private void calculatePopScore() {
        iterCount++;
        if (pop == null || pop.size() == 0) return;
        totalScore = 0;
        averageScore = 0;
        bestScore = 0;
        worstScore = Double.MAX_VALUE;
        bestGene = null;
        worstGene = null;
        for (Chromosome chromosome : pop) {
//            System.out.println(totalScore);
            calculateScore(chromosome);
            totalScore += chromosome.getScore();
            if (chromosome.getScore() > bestScore) {
                bestScore = chromosome.getScore();
                bestGene = chromosome;
            }
            if (chromosome.getScore() < worstScore) {
                worstScore = chromosome.getScore();
                worstGene = chromosome;
            }
            averageScore = totalScore / POP_SIZE;
            averageScore = averageScore > bestScore ? bestScore : averageScore;
        }
        // 将待绘制的折线图数据信息存入 map 中
        bestScoreDataMap.put((double) iterCount, bestScore);
        worstScoreDataMap.put((double) iterCount, worstScore);
        totalScoreDataMap.put((double) iterCount, totalScore);
    }

    // private void calculatePopScore() {
    //     iterCount++;
    //     if (pop==null || pop.size()==0) return;
    //     totalScore=0;
    //     averageScore=0;
    //     bestScore=0;
    //     worstScore=Double.MAX_VALUE;
    //     bestGene=null;
    //     worstGene=null;
    //     for (Chromosome chromosome : pop) {
    //         calculateScore(chromosome);
    //         //System.out.println(chromosome);
    //         totalScore+=chromosome.getScore();
    //         if (chromosome.getScore()>bestScore) {
    //             bestScore=chromosome.getScore();
    //             bestGene= chromosome.getGene();
    //         }
    //         if (chromosome.getScore()<worstScore){
    //             worstScore= chromosome.getScore();
    //             worstGene= chromosome.getGene();
    //         }
    //         averageScore=totalScore/POP_SIZE;
    //         averageScore=averageScore>bestScore?bestScore:averageScore;
    //     }
    //     // 将待绘制的折线图数据信息存入 map 中
    //     bestScoreDataMap.put((double) iterCount,bestScore);
    //     worstScoreDataMap.put((double) iterCount,worstScore);
    //     totalScoreDataMap.put((double) iterCount,totalScore);
    // }


    // TODO: 使用强硬的适应度计算策略 <- 1、在初始化种群的时候 2、在父代交叉的时候
    /**
     * 计算个体适应度 
     * @param chromosome
     */
     private void calculateScore(Chromosome chromosome) {
        ConstraintChecker constraintChecker = new ConstraintChecker();
        
        if (!constraintChecker.isGoodChromosome(chromosome)) {
            throw new RuntimeException("错误：当前染色体出现错误，无法计算");
        }

        double scoreCount = 0d;

                InterferenceCalculator Icalculator = new InterferenceCalculator(chromosome);
                BandwidthConstraint Bcalculator = new BandwidthConstraint(chromosome);

                double interferences = Icalculator.calculateInterferenceLevel(chromosome.getFrequencies(), chromosome.getBandwidths());
                double bandwidthPenalty = Bcalculator.calculateBandwidthPenalty(chromosome.getBandwidths());


                
                // 计算适应度函数 F = w_interference * F_interference + w_bandwidth * F_bandwidth
                scoreCount = (wInterference * interferences) + (wBandwidth * bandwidthPenalty);
    
        // 统计出现负值分数的个数
        if (scoreCount < 0) {
            underZeroCount++;
        }
    
        // 如果适应度分数为0，要做补偿
        scoreCount += ORIGIN_SCORE;
    
        // 最终策略，防止适应度分数小于0，使轮盘赌失效
        // 如果出现大量分数为负数的情况，绘制出的图形会出现异常
        scoreCount = Math.max(0, scoreCount);
        chromosome.setScore(scoreCount);
    }
    
    

    /**
     * 轮盘赌算法  获取较优的父个体
     * @return
     */
    private Chromosome getParentChromosome() {
        // 轮盘赌选中的部分
        double slice = Math.random() * totalScore;
        //System.out.println("轮盘赌 slice:"+slice);
//        System.out.println("轮盘赌 total:"+totalScore);
        double sum = 0d;
//        int targetCount =0;
        for (Chromosome chromosome : pop) {
//            targetCount++;
            sum+=chromosome.getScore();
            // 轮盘赌选中
            if (sum>slice) {
//                System.out.println("当前命中的目标为："+targetCount);
                return chromosome;
            }
//            System.out.println("轮盘赌选中个体为:\n"+chromosome);
        }
        return pop.get(pop.size()-1);
    }

    /**
     * 交叉产生新的子代
     * 这里要注意，一定要避免迭代次数过多导致的程序阻塞
     *
     */
    private void evolve() {
        List<Chromosome> newPop = new ArrayList<>();
        while (newPop.size()<POP_SIZE) {
            Chromosome p1 = getParentChromosome();
            Chromosome p2 = getParentChromosome();
            // genetic 保证生成的子代一定合法
            List<Chromosome> children = Chromosome.genetic(p1, p2);
            ConstraintChecker constraintChecker = new ConstraintChecker();
            for (Chromosome child : children) {
                if (constraintChecker.isGoodChromosome(child)) {
                    newPop.add(child);
                }
            }
        }
        // 保证新生成的子代长度与原父代长度相等
        while (newPop.size()>POP_SIZE) {
            newPop.remove(r.nextInt(newPop.size()));
        }
        pop.clear();
        pop=newPop;
    }

    /**
     * 种群变异
     * 种群变异为了保证变异后子代还是符合条件的
     * 必须对基因组序列进行两两交换
     * 
     * 保证基因组突变后还符合约束
     */
    private void mutation() {
        for (Chromosome chromosome : pop) {
            if (Math.random()<MUTATION_RATE) {
                chromosome.mutation(MAX_MUTATION_NUM,chromosome);
            }
        }
    }
}

    // /**
    //  * 计算两点之间的路径长度
    //  * @param fromX
    //  * @param fromY
    //  * @param destX
    //  * @param destY
    //  * @return
    //  */
    // @Override
    // public double routeLength(double fromX, double fromY, double destX, double destY) {
    //     double deltaX = fromX - destX;
    //     double deltaY = fromY - destY;
    //     return Math.sqrt(deltaX*deltaX + deltaY*deltaY);
    // }

//     /**
//      * 计算基因长短
//      * 因为路径规划问题中，基因长短和车辆个数和受灾点有关系
//      * 所以需要额外开函数计算
//      * @return
//      */
//     private int calculateGeneSize() {
//         return CAR_NUM+POINT_NUM+1;
//     }

// }
