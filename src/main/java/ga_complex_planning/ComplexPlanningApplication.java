package ga_complex_planning;


/**
 * @ClassName ComplexPlanningApplication
 * @Description 主函数
 * @Author fwt
 * @Date 2024/1/29 9:36 下午
 * @Version 1.0
 **/

public class ComplexPlanningApplication {
    public static void main(String[] args) {
        //Map<String, Point> info = Info.getInfo();
        GeneticAlgorithm ga = new GeneticAlgorithm();
        ga.conductGA();

    }
}
