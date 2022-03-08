#include "../Headers/LDP_mechanisms.h"
#include <algorithm>
void mseranges()
{

    double inp_range = 0.0;
    double irr_range = 0.0;
    double hrr_range = 0.0;

    double ihaarht_range = 0.0;

    uint length = 0;
    float h = 0;
    Timer t;

    vector<int> query(2, 0.0);

    double cnt = 0;
    uint d = 22;
    uint step = pow(2, 0);
    uint D = pow(2, d), x = 0;
    uint population = pow(2.0, 26);

    const float e_eps = 3.0;

    Random_Data rand(population, D, 0.33f);

    //std::string filename = "vary_b_ci_" + std::to_string(d) + "_prefix.csv";
    std::string filename = "vary_b_ci_" + std::to_string(d) + ".csv";
    std::ofstream fileptr;
    fileptr.open(filename);
    fileptr << "length"
            << ",";
    fileptr << "count"
            << ",";
    fileptr << "population"
            << ",";
    fileptr << "base"
            << ",";
    fileptr << "OUE"
            << ",";
    fileptr << "HRR"
            << ",";
    fileptr << "TreeOUE"
            << ",";
    fileptr << "TreeHRR"
            << ",";
    fileptr << "TreeOUECI"
            << ",";
    fileptr << "TreeHRRCI"
            << ",";
    fileptr << "HaarHRR"
            << ",";
    fileptr << "TreeOLH"
            << ",";
    fileptr << "TreeOLHCI"
            << "\n";

    fileptr.close();

    vector<AllMechanisms_min> answer, answer_ci;
    for (uint l = 1; l <= D; l++)
    {
        answer.emplace_back();
        answer_ci.emplace_back();
    }

    vector<uint> blist;

    for (uint base = 2; base < D; base++)
    {
        h = log(D) / log(base);
        if ((ceilf(h) == h) && base <= 32)
        {
            blist.push_back(base);
        }
    }
    std::reverse(blist.begin(), blist.end());

    std::cout << blist;
    std::unique_ptr<HaarTransform> ihaar(new HaarTransform(e_eps, D));
    std::unique_ptr<InputHT> iht(new InputHT(e_eps, D));
    InputRR irr(e_eps, D);
    TrueInput inp(D);

    for (uint base : blist)
    {
        //cout << blist ;
        t.start();

        std::unique_ptr<B_ary_tree_min> b_ary(new B_ary_tree_min(0, D - 1, e_eps, base));
        t.stop();
        std::cout << "tree of base =" << base << " , "
                  << "tree allocated in " << t.secs() << std::endl;

        for (short rpt = 0; rpt < 10; rpt++)
        {
            std::cout << "===================" << std::endl;
            std::cout << D << " , "
                      << "rpt=" << rpt << " , "
                      << "base=" << base << std::endl;

            t.start();
            for (uint p = 0; p < population; p++)
            {
                x = rand.sampledata();
                // std::cout << x << std::endl;
                inp.true_input[x] += 1.0f;

                iht->perturb(x);
                irr.perturb(x);
                b_ary->perturb(x);
                ihaar->perturb(x);
            }
            t.stop();
            std::cout << "perturbed in " << t.secs() << std::endl;
            iht->correction();
            irr.correction();
            ihaar->correction();
            inp.normalize();

            t.start();
            b_ary->correction();
            t.stop();
            std::cout << "tree corrected in " << t.secs() << std::endl;
            for (auto it = answer.begin(); it != answer.end(); it++)
                it->clear();

            for (auto it = answer_ci.begin(); it != answer_ci.end(); it++)
                it->clear();

            t.start();
            for (uint i = 0; i < D; i += step)
            {
                query[0] = i;
                inp_range = inp.true_input[i];
                hrr_range = iht->pert_iht[i];
                irr_range = irr.pert_irr[i];

                for (uint j = i + 1; j < D; j++)
                {
                    length = j - i;

                    query[1] = j;

                    inp_range += inp.true_input[j];
                    hrr_range += iht->pert_iht[j];
                    irr_range += irr.pert_irr[j];

                    b_ary->eval_range_itr(query);

                    ihaarht_range = ihaar->eval_range_itr(query);
                    //std:: cout << length << " , " << ihaarht_range << " , " << inp_range << " , " << abs(inp_range-ihaarht_range)  << " , " << irr_range << "\n";
                    //lsanswer[length].length += b_ary->answer.length;
                    answer[length].cnt += 1.0;

                    answer[length].InputRR += pow(irr_range - inp_range, 2);
                    answer[length].HRR += pow(hrr_range - inp_range, 2);
                    answer[length].HaarHRR += pow(ihaarht_range - inp_range, 2);

                    answer[length].TreeRR += pow(b_ary->answer.TreeRR - inp_range, 2);
                    answer[length].TreeHT += pow(b_ary->answer.TreeHT - inp_range, 2);
                    answer[length].TreeOLH += pow(b_ary->answer.TreeOLH - inp_range, 2);
                }
            }
            t.stop();
            std::cout << "queries evaluated in " << t.secs() << std::endl;

            t.start();
            b_ary->constrained_Inference();
            t.stop();
            std::cout << "applied CI " << t.secs() << std::endl;

            t.start();

            for (uint i = 0; i < D; i += step)
            {
                query[0] = i;
                inp_range = inp.true_input[i];

                for (uint j = i + 1; j < D; j++)
                {
                    length = j - i;
                    query[1] = j;

                    inp_range += inp.true_input[j];

                    b_ary->eval_range_itr(query);
                    answer_ci[length].TreeRR += pow(b_ary->answer.TreeRR - inp_range, 2);
                    answer_ci[length].TreeHT += pow(b_ary->answer.TreeHT - inp_range, 2);
                    answer_ci[length].TreeOLH += pow(b_ary->answer.TreeOLH - inp_range, 2);
                }
            }
            t.stop();
            std::cout << "queries after applying CI evaluated in " << t.secs() << std::endl;

            fileptr.open(filename, std::ios::app);
            for (uint i = 1; i < answer.size(); i++)
            {
                cnt = answer[i].cnt;
                if (cnt == 0.0)
                    continue;
                //   fileptr << (answer[i].length / cnt) << ",";
                fileptr << i << ",";
                fileptr << cnt << ",";
                //         fileptr << e_eps << ",";
                fileptr << population << ",";

                fileptr << base << ",";

                fileptr << sqrt(answer[i].InputRR / cnt) << ",";
                fileptr << sqrt(answer[i].HRR / cnt) << ",";

                fileptr << sqrt(answer[i].TreeRR / cnt) << ",";
                fileptr << sqrt(answer[i].TreeHT / cnt) << ",";

                fileptr << sqrt(answer_ci[i].TreeRR / cnt) << ",";
                fileptr << sqrt(answer_ci[i].TreeHT / cnt) << ",";
                fileptr << sqrt(answer[i].HaarHRR / cnt) << ",";
                fileptr << sqrt(answer[i].TreeOLH / cnt) << ",";
                fileptr << sqrt(answer_ci[i].TreeOLH / cnt) << "\n";
            }

            fileptr.close();
            inp.reset_object();
            iht->reset_object();
            irr.reset_object();
            b_ary->reset_object();
            ihaar->reset_object();
        }
    }
}

int main (){

    mseranges();
}