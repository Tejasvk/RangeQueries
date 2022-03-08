#include "../Headers/LDP_mechanisms.h"
#include <algorithm>

/*
Use this function when you want to find length wise MSE for various epsilon values.
*/
void avg_error_eps()
{

/// Setting this macro to 1 includes OLH method which is compute intensive.
/// We do not include it by default.
#define full 0
    double inp_range = 0.0;
    double irr_range = 0.0;
    double hrr_range = 0.0;
    double ihaarht_range = 0.0;

    uint length = 0;
    float h = 0;
    Timer t;

    uint d = 20;
    // Use a step size other than 1 when D is very large say > 2**20 in order to skip some queries.
    // Otherwise, let step size be 1.
    uint step = pow(2, 15);
    uint D = pow(2, d), x = 0;
    uint population = pow(2.0, 26);
    std::vector<int> query(2, 0);

    Random_Data rand(population, D, 0.4f);

    std::string filename = "avg_case_error_" + std::to_string(d) + ".csv";
    //std::string filename = "avg_case_error_" + std::to_string(d) + "_pre.csv";
    std::ofstream fileptr;

    fileptr.open(filename);

    fileptr << "e_eps"
            << ",";
    fileptr << "population"
            << ",";
    fileptr << "base"
            << ",";
    fileptr << "OUE"
            << ",";
    fileptr << "HRR"
            << ",";
    fileptr << "TreeOUECI"
            << ",";
    fileptr << "TreeHRRCI"
            << ",";
    fileptr << "HaarHRR"
            << ",";
#if full == 1
    fileptr << "TreeOLHCI";
#endif
    fileptr << "\n";

    fileptr.close();

    AllMechanisms_min avg_ans_ci;

    vector<uint> blist;
    // We consider base values up to 256.
    for (uint base = 2; base < D; base++)
    {
        h = log(D) / log(base);
        if (ceilf(h) == h && base <= 256)
        {
            blist.push_back(base);
        }
    }
    std::reverse(blist.begin(), blist.end());

    std::cout << blist;
    TrueInput inp(D);
    vector<float> eps_list = {0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.1f, 1.2f, 1.4f};

    std::reverse(eps_list.begin(), eps_list.end());
    for (float eps : eps_list)
    {
        float e_eps = power(eps);
        std::unique_ptr<HaarTransform> ihaar(new HaarTransform(e_eps, D));
        std::unique_ptr<InputHT> iht(new InputHT(e_eps, D));
        InputRR irr(e_eps, D);

        for (uint base : blist)
        {
            //cout << blist ;
            t.start();

            std::unique_ptr<B_ary_tree_min> b_ary(new B_ary_tree_min(0, D - 1, e_eps, base));
            t.stop();
            std::cout << "tree of base =" << base << " , "
                      << "tree allocated in " << t.secs() << std::endl;

            for (short rpt = 0; rpt < 5; rpt++)
            {
                std::cout << "===================" << std::endl;
                std::cout << D << " , "
                          << "eps=" << eps << " , "
                          << "rpt=" << rpt << " , "
                          << "base=" << base << std::endl;

                t.start();
                // Perturbation simulation.
                for (uint p = 0; p < population; p++)
                {
                    x = rand.sampledata();
                    inp.true_input[x] += 1.0;
                    iht->perturb(x);
                    irr.perturb(x);
                    b_ary->perturb(x);
                    ihaar->perturb(x);
                }
                t.stop();
                std::cout << "perturbed in " << t.secs() << std::endl;

                iht->correction();
                irr.correction3();
                ihaar->correction();
                inp.normalize();

                t.start();
                b_ary->correction();
                t.stop();
                std::cout << "tree corrected in " << t.secs() << std::endl;

                t.start();
                b_ary->constrained_Inference();
                t.stop();
                std::cout << "applied CI " << t.secs() << std::endl;

                t.start();

                for (uint i = 0; i < D; i += step)
                {
                    query[0] = i;
                    inp_range = inp.true_input[i];
                    hrr_range = iht->pert_iht[i];
                    irr_range = irr.pert_irr[i];

                    for (uint j = i + 1; j < D; j++)
                    {
                        query[1] = j;

                        inp_range += inp.true_input[j];

                        hrr_range += iht->pert_iht[j];
                        avg_ans_ci.HRR += pow(hrr_range - inp_range, 2);
                        irr_range += irr.pert_irr[j];
                        avg_ans_ci.InputRR += pow(irr_range - inp_range, 2);

                        ihaarht_range = ihaar->eval_range_itr(query);
                        avg_ans_ci.HaarHRR += pow(ihaarht_range - inp_range, 2);

                        b_ary->eval_range_itr(query);
                        avg_ans_ci.TreeRR += pow(b_ary->answer.TreeRR - inp_range, 2);
                        avg_ans_ci.TreeHT += pow(b_ary->answer.TreeHT - inp_range, 2);
#if full == 1
                        avg_ans_ci.TreeOLH += pow(b_ary->answer.TreeOLH - inp_range, 2);
#endif
                        avg_ans_ci.cnt += 1.0;
                    }
                }
                t.stop();
                std::cout << "total queries executed = " << (uint)avg_ans_ci.cnt << std::endl;
                std::cout << "queries after applying CI evaluated in " << t.secs() << std::endl;

                fileptr.open(filename, std::ios::app);

                double cnt = avg_ans_ci.cnt;
                fileptr << e_eps << ",";
                fileptr << population << ",";
                fileptr << base << ",";
                fileptr << sqrt(avg_ans_ci.InputRR / cnt) << ",";
                fileptr << sqrt(avg_ans_ci.HRR / cnt) << ",";

                fileptr << sqrt(avg_ans_ci.TreeRR / cnt) << ",";
                fileptr << sqrt(avg_ans_ci.TreeHT / cnt) << ",";
                fileptr << sqrt(avg_ans_ci.HaarHRR / cnt);
#if full == 1
                fileptr << "," << sqrt(avg_ans_ci.TreeOLH / cnt);
#endif
                fileptr << "\n";
                fileptr.close();
                avg_ans_ci.clear();

                inp.reset_object();
                irr.reset_object();
                iht->reset_object();
                b_ary->reset_object();
                ihaar->reset_object();
            }
        }
    }
}

/*
Use this function when you want to find length wise MSE for various distribution shapes.
*/
void avg_error_dbn()
{

#define full 0

    double inp_range = 0.0;
    double irr_range = 0.0;
    double hrr_range = 0.0;
    double ihaarht_range = 0.0;

    float e_eps = 3.0;
    uint length = 0;
    float h = 0;
    Timer t;
    uint d = 16; /// Change the D size.
    uint step = pow(2, 11);
    uint D = pow(2, d), x = 0;
    uint base = 4;
    uint population = pow(2.0, 26);
    std::vector<int> query(2, 0.0);

    //std::string filename = "avg_case_error_" + std::to_string(d) + ".csv";
    std::string filename = "avg_case_error_" + std::to_string(d) + "_dbn.csv";

    std::ofstream fileptr;
    fileptr.open(filename);
    fileptr << "prob"
            << ",";
    fileptr << "e_eps"
            << ",";
    fileptr << "population"
            << ",";
    fileptr << "base"
            << ",";
    fileptr << "OUE"
            << ",";
    fileptr << "HRR"
            << ",";
    fileptr << "TreeOUECI"
            << ",";
    fileptr << "TreeHRRCI"
            << ",";
    fileptr << "HaarHRR";
#if full == 1
    fileptr << ",";
    fileptr << "TreeOLH"
            << ",";
    fileptr << "TreeOLHCI"
#endif
        fileptr
            << "\n";

    fileptr.close();
    AllMechanisms_min avg_ans_ci;

    TrueInput inp(D);
    InputRR irr(e_eps, D);

    std::unique_ptr<InputHT> iht(new InputHT(e_eps, D));

    std::unique_ptr<HaarTransform> ihaar(new HaarTransform(e_eps, D));

    vector<float> prob_list = {0.0f, 0.6f, 0.2f, 0.3f, 0.4f, 0.5f, 0.1f, 0.7f, 0.8f, 0.9f, 1.0f};
    for (float prob : prob_list)
    {
        Random_Data rand(population, D, prob);

        t.start();

        std::unique_ptr<B_ary_tree_min> b_ary(new B_ary_tree_min(0, D - 1, e_eps, base));
        t.stop();
        std::cout << "tree of base =" << base << " , "
                  << "tree allocated in " << t.secs() << std::endl;

        for (short rpt = 0; rpt < 5; rpt++)
        {
            std::cout << "===================" << std::endl;
            std::cout << D << " , "
                      << "rpt=" << rpt << " , "
                      << "base=" << base << std::endl;

            t.start();
            for (uint p = 0; p < population; p++)
            {
                x = rand.sampledata();
                inp.true_input[x] += 1.0;
                iht->perturb(x);
                irr.perturb(x);
                b_ary->perturb(x);
                ihaar->perturb(x);
            }
            t.stop();
            std::cout << "perturbed in " << t.secs() << std::endl;
            iht->correction();
            irr.correction3();
            ihaar->correction();
            inp.normalize();

            t.start();
            b_ary->correction();
            t.stop();
            std::cout << "tree corrected in " << t.secs() << std::endl;
            t.start();
            b_ary->constrained_Inference();
            t.stop();
            std::cout << "applied CI " << t.secs() << std::endl;

            t.start();

            for (uint i = 0; i < D; i += step)
            {
                query[0] = i;
                inp_range = inp.true_input[i];
                hrr_range = iht->pert_iht[i];
                irr_range = irr.pert_irr[i];

                for (uint j = i + 1; j < D; j++)
                {

                    inp_range += inp.true_input[j];
                    hrr_range += iht->pert_iht[j];
                    irr_range += irr.pert_irr[j];

                    query[1] = j;
    
                    avg_ans_ci.InputRR += pow(irr_range - inp_range, 2);
                    avg_ans_ci.HRR += pow(hrr_range - inp_range, 2);
                    ihaarht_range = ihaar->eval_range_itr(query);
                    avg_ans_ci.HaarHRR += pow(ihaarht_range - inp_range, 2);

                    b_ary->eval_range_itr(query);
                    avg_ans_ci.TreeRR += pow(b_ary->answer.TreeRR - inp_range, 2);
                    avg_ans_ci.TreeHT += pow(b_ary->answer.TreeHT - inp_range, 2);
#if full == 1

                    avg_ans_ci.TreeOLH += pow(b_ary->answer.TreeOLH - inp_range, 2);
#endif
                    avg_ans_ci.cnt += 1.0;
                }
            }
        }
        t.stop();
        std::cout << "queries executed " << (uint)avg_ans_ci.cnt << std::endl;
        std::cout << "queries after applying CI evaluated in " << t.secs() << std::endl;

        fileptr.open(filename, std::ios::app);

        double cnt = avg_ans_ci.cnt;

        fileptr << prob << ",";
        fileptr << e_eps << ",";
        fileptr << population << ",";
        fileptr << base << ",";
        fileptr << sqrt(avg_ans_ci.InputRR / cnt) << ",";
        fileptr << sqrt(avg_ans_ci.HRR / cnt) << ",";
        fileptr << sqrt(avg_ans_ci.TreeRR / cnt) << ",";
        fileptr << sqrt(avg_ans_ci.TreeHT / cnt) << ",";
        fileptr << sqrt(avg_ans_ci.HaarHRR / cnt);
#if full == 1
        fileptr << "," << sqrt(avg_ans_ci.TreeOLH / cnt);
#endif
        fileptr << "\n";
        fileptr.close();

        avg_ans_ci.clear();
        inp.reset_object();

        irr.reset_object();
        iht->reset_object();
        b_ary->reset_object();
        ihaar->reset_object();
    }
}


int main()
{

    //avg_error_eps();
    //avg_error_dbn();
}