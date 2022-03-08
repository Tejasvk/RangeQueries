#include "../Headers/LDP_mechanisms.h"

void quantile_dbn()
{
#define raw 1
    uint population = pow(2.0, 20);
    float e_eps = 3.0f;
    float h = 16;
    double inp_range = 0.0;
    double irr_range = 0.0;
    double ihaar_range = 0.0;  
    uint x = 0;
    uint D = pow(2, 10);

    std::string filename = "quantile_prob_" + std::to_string(int(D)) + ".csv";
    std::ofstream fileptr;
    fileptr.open(filename);

    fileptr << "prob"
            << ",";
    fileptr << "quantile"
            << ",";
    fileptr << "population"
            << ",";
    fileptr << "e_eps"
            << ",";
    fileptr << "base"
            << ",";
#if raw == 1
    fileptr << "True"
            << ",";
#endif
    fileptr << "OUE"
            << ",";
    fileptr << "TreeOUECI"
            << ",";
    fileptr << "TreeHRRCI"
            << ",";
    fileptr << "HaarHRR";
#if raw == 1
    fileptr << ",";
    fileptr << "QTrue"
            << ",";
    fileptr << "QOUE"
            << ",";
    fileptr << "QTreeOUECI"
            << ",";
    fileptr << "QTreeHRRCI"
            << ",";
    fileptr << "QHaarHRR";
#endif
    fileptr << "\n";

    fileptr.close();

    vector<int> query(2, 0);
    uint base = 4;
    double quant = 0.1;
    std::unique_ptr<B_ary_tree_min> b_ary(new B_ary_tree_min(0, D - 1, e_eps, base));

    std::unique_ptr<HaarTransform> ihaar(new HaarTransform(e_eps, D));

    float prob = 0.1;

    while (prob < 1.0)
    {
        Random_Data rand(population, D, prob);
        TrueInput inp(D);
        InputRR irr(e_eps, D);

        for (short rpt = 0; rpt < 5; rpt++)
        {

            for (uint p = 0; p < population; p++)
            {
                x = rand.sampledata();
                inp.true_input[x] += 1.0f;
                irr.perturb(x);
                b_ary->perturb(x);
                ihaar->perturb(x);
            }

            inp.normalize();
            irr.correction3();
            b_ary->correction();
            b_ary->constrained_Inference();
            ihaar->correction();
            quant = 0.1;
            while (quant < 1.0)
            {

                std::vector<double> ans(5, 0.0);
                std::vector<double> qindex(5, 0.0);

                bool flags[5] = {true, true, true, true, true};

                inp_range = inp.true_input[0];
                irr_range = irr.pert_irr[0];

                uint len_hh = 0, len_haar = 0;
                for (uint j = 1; j < D; j++)
                {
                    if (flags[0])
                    {
                        inp_range += inp.true_input[j];
                        if (inp_range >= quant)
                        {
                            ans[0] = inp_range;
                            flags[0] = false;
                            qindex[0] = j;
                        }
                    }

                    if (flags[1])
                    {
                        irr_range += irr.pert_irr[j];
                        if (irr_range >= quant)
                        {
                            ans[1] = irr_range;
                            flags[1] = false;
                            qindex[1] = j;
                        }
                    }
                    if (flags[2])
                    {
                        query[0] = 0;
                        query[1] = j;

                        b_ary->eval_range_itr(query);
                        if (b_ary->answer.TreeRR >= quant)
                        {
                            ans[2] = b_ary->answer.TreeRR;
                            flags[2] = false;
                            len_hh = b_ary->answer.length;
                            qindex[2] = j;
                        }
                    }
                    if (flags[3])
                    {
                        query[0] = 0;
                        query[1] = j;

                        b_ary->eval_range_itr(query);
                        if (b_ary->answer.TreeHT >= quant)
                        {
                            ans[3] = b_ary->answer.TreeHT;
                            flags[3] = false;
                            qindex[3] = j;
                        }
                    }
                    if (flags[4])
                    {
                        query[0] = 0;
                        query[1] = j;

                        ihaar_range = ihaar->eval_range_itr(query);
                        if (ihaar_range >= quant)
                        {
                            // std::cout << ihaar_range <<  "," << quant << "\n";
                            ans[4] = ihaar_range;
                            flags[4] = false;
                            len_haar = ihaar->no_nodes;
                            qindex[4] = j;
                        }
                    }

                    if (!(flags[0] || flags[1] || flags[2] || flags[3] || flags[4]))
                        break;
                }

                //std::cout << "prob=" << prob << " ," << true_ind << " , " << 100.0 * abs((float)true_ind - flat_ind) / true_ind << " , " << 100.0 * abs((float)hh_ind - true_ind) / true_ind << " , " << 100.0 * abs((float)haar_ind - true_ind) / true_ind << "\n";

                if ((ans[0] != 0.0) && (ans[1] != 0.0) && (ans[2] != 0.0) && (ans[3] != 0.0) && (ans[4] != 0.0))
                {
                    fileptr.open(filename, std::ios::app);
                    fileptr << prob << ",";
                    fileptr << quant << ",";
                    fileptr << population << ",";
                    fileptr << e_eps << ",";
                    fileptr << base << ",";

#if raw == 0
                    fileptr << abs(qindex[0] - qindex[1]) << ",";
                    fileptr << abs(qindex[0] - qindex[2]) << ",";
                    fileptr << abs(qindex[0] - qindex[3]) << ",";
                    fileptr << abs(qindex[0] - qindex[4]) << "\n";
#else
                    fileptr << qindex[0] << ",";
                    fileptr << qindex[1] << ",";
                    fileptr << qindex[2] << ",";
                    fileptr << qindex[3] << ",";
                    fileptr << qindex[4] << ",";

                    fileptr << ans[0] << ",";
                    fileptr << ans[1] << ",";
                    fileptr << ans[2] << ",";
                    fileptr << ans[3] << ",";
                    fileptr << ans[4] << "\n";
#endif

                    fileptr.close();
                }
                quant += 0.1f;
            }
#if raw == 0
            fileptr.open(filename, std::ios::app);
            fileptr << "\n";
            fileptr.close();
#endif
            b_ary->reset_object();
            ihaar->reset_object();
            irr.reset_object();
            inp.reset_object();
        }
        prob += 0.1;
    }
}

int main()
{

    quantile_dbn();

}