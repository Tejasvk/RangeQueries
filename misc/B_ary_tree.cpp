/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "../Headers/B_ary_tree.h"

void B_ary_tree::reset_object()
{

    for (size_t i = 0; i < uirrarray.size(); i++)
    {
        uirrarray[i].reset();
        olharray[i].reset();
        lhtarray[i].reset();
    }
    uirrarray.erase(uirrarray.begin(), uirrarray.end());
    olharray.erase(olharray.begin(), olharray.end());
    lhtarray.erase(lhtarray.begin(), lhtarray.end());
 
    level = 1;

    //inphaar= std::make_unique<HaarTransform>(e_eps, r_range+1);

    while (level <= height)
    {
        uirrarray.push_back(std::make_unique<InputRR>(e_eps, (uint)pow(b, level)));
        olharray.push_back(std::make_unique<InputOLH>(e_eps, (uint)pow(b, level)));
        lhtarray.push_back(std::make_unique<InputHT>(e_eps, (uint)pow(b, level)));
 
        level++;
    }

    for (auto it = arraytree.begin(); it != arraytree.end(); it++)
    {
        (*it)->clearcounts();
    }
}

B_ary_tree::B_ary_tree(int l, int r, float e_eps, int b) : l_range(l), r_range(r), e_eps(e_eps), b(b)
{

    height = ceil(log(r) / log(b));
    std::cout << "height="
              << "," << height << std::endl;
    std::cout << r << std::endl;

    buildArrayTree();
    std::cout << "Tree allocated!" << std::endl;

    level = 1;
    //vector<int> levelnodesids;
    //inphaar= std::make_unique<HaarTransform>(e_eps, r_range+1);
    std::cout << height << std::endl;
    while (level <= height)
    {
        //cout << "allocating level objects " << level << "," << (uint) pow(b, level) << endl;
        uirrarray.push_back(std::make_unique<InputRR>(e_eps, (uint)pow(b, level)));
        olharray.push_back(std::make_unique<InputOLH>(e_eps, (uint)pow(b, level)));
        lhtarray.push_back(std::make_unique<InputHT>(e_eps, (uint)pow(b, level)));
        level++;
    }
    std::cout << "allocated all level objects" << std::endl;
    uni_level_sampler = std::uniform_int_distribution<int>(1, height);
    //uni_level_sampler_dht = std::uniform_int_distribution<int>(1, height);
    buildIndexCache3();
    std::cout << "Level cache allocated!" << std::endl;
    //cout << endl << "height=" << height << endl;
    //generator.seed(time(NULL));
}

void B_ary_tree::buildIndexCache3()
{

    int j = 0;
    int l = 0;
    int r = 0;
    trueindexcache.reserve(r_range);

    for (level = 1; level <= height; level++)
    {
        j = 0;
        l = 0;
        while (j < level)
            l += pow(b, j++);
        r = l + pow(b, level);
        j = 0;
        //cout << "-----------------" << endl;
        while (l < r)
        {
            arraytree[l]->hid = j++;
            //cout << arraytree[l]->hid <<"," <<arraytree[l]->range;
            l++;
        }
        //        cout << level << endl;
    }
    for (int x = 0; x <= r_range; x++)
    {
        n = arraytree[leafoffset + x];
        vector<uint> path;
        path.reserve(height);
        while (n->parent)
        {
            path.emplace_back(n->hid);
            n = n->parent;
        }
        //path.pop_back();
        //cout << path;
        //std::reverse(path.begin(),path.end());
        trueindexcache.emplace_back(path);
    }
}

vector<node_b *> B_ary_tree::buildArrayTree()
{
    std::queue<std::tuple<std::vector<float>, node_b *>> q;
    float l_temp = l_range;
    float r_temp = r_range;

    float l_temp_lvl;
    float r_temp_lvl;
    int level = 0;

    phi = {l_temp, r_temp};
    int b_split = 0, b_split_prev = 0;
    int id = 0, i = 0;
    node_b *root = NULL;
    n = new node_b(l_temp, r_temp, id, b);
    //cout << "size " << sizeof (*n) << endl;
    n->level = 0;
    q.push(std::make_tuple(phi, n));
    arraytree.emplace_back(n);
    while (!q.empty())
    {
        std::tie(phi, root) = q.front();
        q.pop();

        if (phi[0] > phi[1])
            continue;

        l_temp = phi[0];
        r_temp = phi[1];
        //cout << root->id << endl;
        b_split = ceil((r_temp - l_temp) / b);
        if (b_split != b_split_prev)
            level++;
        //cout << b_split << "," << level <<  endl;
        l_temp_lvl = l_temp;
        r_temp_lvl = l_temp_lvl + b_split - 1.0;
        i = 0;
        //cout << "----------" << endl;
        while ((i < b) && (l_temp_lvl <= r_temp_lvl))
        {
            n = new node_b(l_temp_lvl, r_temp_lvl, ++id, b);
            n->parent = root;
            n->level = level;
            root->children[i] = n;

            //q.push(std::make_tuple(vector<float>{l_temp_lvl, r_temp_lvl}, n));
            q.emplace(vector<float>{l_temp_lvl, r_temp_lvl}, n);
            arraytree.emplace_back(n);

            l_temp_lvl = r_temp_lvl + 1.0;
            r_temp_lvl = l_temp_lvl + b_split - 1.0;
            i++;
        }
        b_split_prev = b_split;
    }
    leafoffset = arraytree.size() - r_range - 1;
    //cout << endl << q.size() << endl;
    std::queue<std::tuple<std::vector<float>, node_b *>> empty;
    std::swap(q, empty);
    //cout << "=============" << endl;
    return arraytree;
}

void B_ary_tree::perturb(const uint &x)
{

    level = uni_level_sampler(generator);
    trueindex = trueindexcache[x][height - level];
    olharray[level - 1]->perturb(trueindex);
    lhtarray[level - 1]->perturb(trueindex);
    // lhtarray[level - 1]->perturb_more_coefs(trueindex);
    uirrarray[level - 1]->perturb(trueindex);
}

void B_ary_tree::correction()
{

    arraytree[0]->pert_freq_irr = 1.0;
    arraytree[0]->level_olh = 1.0;
    arraytree[0]->level_ht = 1.0;
    //arraytree[0]->level_ht_ps = 1.0;

    level = 1;
    //cout << "====================" << endl;
    //cout << "correction" << endl;

    int j = 0;
    int l = 0;
    int r = 0;
    // cout << " correction started!" << endl;
    for (level = 1; level <= height; level++)
    {
        //cout << levelnodes;

        uirrarray[level - 1]->correction();
        olharray[level - 1]->correction();
        lhtarray[level - 1]->correction();
        //cout << "====================" << endl;
        //  cout << lhaararray[level]->pert_dht_ps;
        //c_sum.clear();
        j = 0;
        l = 0;
        while (j < level)
            l += pow(b, j++);
        r = l + pow(b, level);
        //cout << "-----------------" << endl;
        j = 0;
        //cout << l << "," << r << "," << level << endl;
        while (l < r)
        {

            arraytree[l]->pert_freq_irr = uirrarray[level - 1]->pert_irr[j];
            arraytree[l]->level_olh = olharray[level - 1]->pert_olh[j];
            arraytree[l]->level_ht = lhtarray[level - 1]->pert_iht[j];
            //arraytree[l]->level_ht_ps = lhtarray[level - 1]->pert_iht_multi_coef[j];

            j++;
            l++;
        }
    }
}

AllMechanisms B_ary_tree::eval_range_itr(vector<int> &range)
{

    if (arraytree[0] == NULL)
        return answer;
    answer.clear();
    vector<int> phi_temp(2, 0);
    phi_temp[0] = l_range;
    phi_temp[1] = r_range;
    //vector<int>phi2 = phi;

    n = arraytree[0];
    nodeStack.emplace(phi_temp, n);
    //cout << "==============" << endl;
    //cout <<"query" << range;
    //cout << "--------------" << endl;
    while (!nodeStack.empty())
    {
        std::tie(phi_temp, n) = nodeStack.top();
        nodeStack.pop();
        if ((phi_temp[0] > phi_temp[1]) || (!n))
            continue;

        phi_temp[0] = std::max(n->range[0], range[0]);
        phi_temp[1] = std::min(n->range[1], range[1]);

        if ((phi_temp[0] == n->range[0]) && (phi_temp[1] == n->range[1]))
        {
            answer.TreeRR += n->pert_freq_irr;
            answer.TreeHT += n->level_ht;
            //  answer.level_ht_ps += n->level_ht_ps;
            answer.TreeOLH += n->level_olh;
            answer.length += 1;
            //cout << range;
            continue;
        }
        for (int i = 0; i < b; i++)
        {
            if (!(n->children[i]))
                break;
            //phi2 = {max(n->children[i]->range[0], phi[0]), min(n->children[i]->range[1], phi[1])};
            //if (phi2[0] <= phi2[1]) {
            if ((std::max(n->children[i]->range[0], phi_temp[0])) <= std::min(n->children[i]->range[1], phi_temp[1]))
            {
                nodeStack.emplace(phi_temp, n->children[i]);
                //cout << n->children[i]->range;
            }
        }
    }

    return answer;
}

void B_ary_tree::constrained_Inference()
{
    weighted_averaging(arraytree[0]);
    mean_Consistency(arraytree[0]);
    //normalize_Tree();
}

void B_ary_tree::weighted_averaging(node_b *root)
{
    if (!root)
        return;
    for (int i = 0; i < b; i++)
        weighted_averaging(root->children[i]);

    if ((height - root->level) > 0)
    {
        float root_weight = (pow(b, height - root->level) - pow(b, height - root->level - 1.0)) / (pow(b, height - root->level) - 1.0);
        float children_weight = 1.0 - root_weight;
        //cout << root_weight  << "," << root->level << " , " << root->range;
        c_sum.clear();
        for (int i = 0; i < b; i++)
        {
            if (root->children[i])
            {
                c_sum.TreeRR += root->children[i]->pert_freq_irr;
                c_sum.TreeHT += root->children[i]->level_ht;
                c_sum.TreeOLH += root->children[i]->level_olh;
            }
        }

        root->pert_freq_irr = (root->pert_freq_irr * root_weight) + (c_sum.TreeRR * children_weight);
        root->level_ht = (root->level_ht * root_weight) + (c_sum.TreeHT * children_weight);
        root->level_olh = (root->level_olh * root_weight) + (c_sum.TreeOLH * children_weight);
    }
}

void B_ary_tree::mean_Consistency(node_b *root)
{
    if (!root)
        return;

    c_sum.clear();
    for (int i = 0; i < b; i++)
    {
        if (root->children[i])
        {
            c_sum.length += 1.0;
            c_sum.TreeRR += root->children[i]->pert_freq_irr;
            c_sum.TreeHT += root->children[i]->level_ht;
            c_sum.TreeOLH += root->children[i]->level_olh;
        }
    }

    for (int i = 0; i < b; i++)
    {
        if (root->children[i])
        {
            //cout << abs(-c_sum.sm_rr_u + root->pert_freq_irr_u) / c_sum.length << endl;

            root->children[i]->pert_freq_irr += ((root->pert_freq_irr - c_sum.TreeRR) / c_sum.length);
            root->children[i]->level_ht += ((root->level_ht - c_sum.TreeHT) / c_sum.length);
            root->children[i]->level_olh += ((root->level_olh - c_sum.TreeOLH) / c_sum.length);
        }
    }
    for (int i = 0; i < b; i++)
    {
        if (root->children[i])
        {
            mean_Consistency(root->children[i]);
        }
    }
}
