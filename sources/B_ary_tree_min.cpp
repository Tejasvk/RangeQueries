

#include "../Headers/B_ary_tree_min.h"


/*
We reset counters for all the objects but retains the caches.
Use this function when you want to reuse the object after clearing the counters.
*/
void B_ary_tree_min::reset_object()
{

    uirrarray.clear();
    uirrarray.erase(uirrarray.begin(), uirrarray.end());

#if OLH == 1
    olharray.clear();
    olharray.erase(olharray.begin(), olharray.end());
#endif

    lhtarray.clear();
    lhtarray.erase(lhtarray.begin(), lhtarray.end());
    
    // reallocate the HRR and OUE level objects.
    level = 1;
    while (level <= height)
    {
        uirrarray.push_back(std::make_unique<InputRR>(e_eps, (uint)pow(b, level)));
        lhtarray.push_back(std::make_unique<InputHT>(e_eps, (uint)pow(b, level)));

#if OLH == 1
        olharray.push_back(std::make_unique<InputOLH>(e_eps, (uint)pow(b, level)));
#endif
        level++;
    }

 // Clear tree level counts.
    for (auto it = arraytree.begin(); it != arraytree.end(); it++)
    {
        (*it)->clearcounts();
    }

    std::cout << "tree object reset" << std::endl;
}

B_ary_tree_min::B_ary_tree_min(int l, int r, float e_eps, int b) : l_range(l), r_range(r), e_eps(e_eps), b(b)
{
    std::cout << "inside tree constructor" << std::endl;
    height = ceil(log(r) / log(b));
    //cout << "height=" << "," << height << endl;
    buildArrayTree();
    level = 1;
    //vector<int> levelnodesids;
    //cout << "height of the tree is =" <<height << endl;
    // allocating the tree level objects.
    while (level <= height)
    {
        //cout << "allocating level objects " << level << "," << (uint) pow(b, level) << endl;

        uirrarray.push_back(std::make_unique<InputRR>(e_eps, (uint)pow(b, level)));
        lhtarray.push_back(std::make_unique<InputHT>(e_eps, (uint)pow(b, level)));

#if OLH == 1
        olharray.push_back(std::make_unique<InputOLH>(e_eps, (uint)pow(b, level)));
        std:: cout << "adding OLH" << std::endl;
#endif
        level++;
    }
    //cout << "allocated all level objects" << endl;

    uni_level_sampler = std::uniform_int_distribution<int>(1, height);
    buildIndexCache();
    //  std::cout << "Level cache allocated!" << std::endl;
}

/*
The index cache method facililates retrieving the signal index for a given item & level pair in constant time.
*/
void B_ary_tree_min::buildIndexCache()
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
        while (l < r)
        {
            arraytree[l]->hid = j++;
            l++;
        }
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
        trueindexcache.emplace_back(path);
    }
}
/*
We build a B-adic tree. The root is stored at  arraytree[0]. 
Note that, though the tree is stored in a vector, we also maintain parent, children links for convenience.  
*/
void B_ary_tree_min::buildArrayTree()
{
    std::queue<std::tuple<std::vector<float>, node_b_l *>> q;
    float l_temp = l_range;
    float r_temp = r_range;

    float l_temp_lvl;
    float r_temp_lvl;
    int level = 0;

    phi = {l_temp, r_temp};
    int b_split = 0, b_split_prev = 0;
    int id = 0, i = 0;
    node_b_l *root = NULL;
    n = new node_b_l(l_temp, r_temp, id, b);
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
        b_split = ceil((r_temp - l_temp) / b);
        if (b_split != b_split_prev)
            level++;

        l_temp_lvl = l_temp;
        r_temp_lvl = l_temp_lvl + b_split - 1.0f;
        i = 0;
        while ((i < b) && (l_temp_lvl <= r_temp_lvl))
        {

            n = new node_b_l(l_temp_lvl, r_temp_lvl, ++id, b);
            n->parent = root;
            n->level = level;
            root->children[i] = n;
            q.emplace(vector<float>{l_temp_lvl, r_temp_lvl}, n);
            arraytree.emplace_back(n);

            l_temp_lvl = r_temp_lvl + 1.0;
            r_temp_lvl = l_temp_lvl + b_split - 1.0;
            i++;
        }
        b_split_prev = b_split;
    }
    leafoffset = arraytree.size() - r_range - 1;

    std::queue<std::tuple<std::vector<float>, node_b_l *>> empty;
    std::swap(q, empty);
    std::cout << "tree allocated" << std::endl;

}


/*
Call the function in the main simulation loop. 
*/
void B_ary_tree_min::perturb(const uint &x)
{
    level = uni_level_sampler(generator);

    trueindex = trueindexcache[x][height - level];
    uirrarray[level - 1]->perturb(trueindex);
    lhtarray[level - 1]->perturb(trueindex);
#if OLH == 1
    olharray[level - 1]->perturb(trueindex);
#endif
}

/*
Corrects the tree levelwise.
*/
void B_ary_tree_min::correction()
{

    arraytree[0]->pert_freq_irr = 1.0f;
    arraytree[0]->level_ht = 1.0f;
#if OLH == 1
    arraytree[0]->level_olh = 1.0f;
#endif
    level = 1;

    int j = 0;
    int l = 0;
    int r = 0;

    for (level = 1; level <= height; level++)
    {
        uirrarray[level - 1]->correction();
        lhtarray[level - 1]->correction();
#if OLH == 1
        olharray[level - 1]->correction();
#endif
        j = 0;
        l = 0;
        while (j < level)
            l += pow(b, j++);
        r = l + pow(b, level);
        j = 0;
        while (l < r)
        {
            arraytree[l]->pert_freq_irr = uirrarray[level - 1]->pert_irr[j];
            arraytree[l]->level_ht = lhtarray[level - 1]->pert_iht[j];

#if OLH == 1
            arraytree[l]->level_olh = olharray[level - 1]->pert_olh[j];
#endif
            j++;
            l++;
        }
    }

    //  cout << "corrected" << endl;
}

/*
We walk over the tree and iterate a given query. 
*/
void B_ary_tree_min::eval_range_itr(vector<int> &range)
{

    //if (arraytree[0] == NULL)
    //    return;
    answer.clear();
    phi_range[0] = l_range;
    phi_range[1] = r_range;
    n = arraytree[0];

    nodeStack.emplace(phi_range, n);
    // cout << "==============" << endl;
    //cout <<"query" << range;
    //cout << "--------------" << endl;
    int i = 0;
    while (!nodeStack.empty())
    {
        std::tie(phi_range, n) = nodeStack.top();
        nodeStack.pop();
        if ((phi_range[0] > phi_range[1]) || (!n))
            continue;

        phi_range[0] = std::max(n->range[0], range[0]);
        phi_range[1] = std::min(n->range[1], range[1]);

        if ((phi_range[0] == n->range[0]) && (phi_range[1] == n->range[1]))
        {

            answer.TreeRR += n->pert_freq_irr;
            answer.TreeHT += n->level_ht;
#if OLH == 1
            answer.TreeOLH += n->level_olh;
#endif
            answer.length += 1;
            continue;
        }
        for (i = 0; i < b; i++)
        {
            if (!(n->children[i]))
                break;
            phi2[0] = std::max(n->children[i]->range[0], phi_range[0]);
            phi2[1] = std::min(n->children[i]->range[1], phi_range[1]);
            if (phi2[0] <= phi2[1])
            {
                nodeStack.emplace(phi_range, n->children[i]);
            }
        }
    }

}
/*
This is a wrapper that applies the 2 stage consistency framework. 
*/
void B_ary_tree_min::constrained_Inference()
{

    weighted_averaging_itr(arraytree[0]);
    mean_Consistency_itr(arraytree[0]);
}

/*
These is an iterative version of the mean consistency step.
*/
void B_ary_tree_min::mean_Consistency_itr(node_b_l *root)
{

    int l = 0;
    int r = 0;
    for (level = 0; level < height; level++)
    {

        j = 0;
        l = 0;
        while (j < level)
            l += pow(b, j++);
        r = l + pow(b, level);
        while (l < r)
        {
            c_sum.clear();
            for (int i = 0; i < b; i++)
            {
                if (!arraytree[l]->children[i])
                    break;
                c_sum.length += 1.0;
                c_sum.TreeRR += arraytree[l]->children[i]->pert_freq_irr;
                c_sum.TreeHT += arraytree[l]->children[i]->level_ht;
#if OLH == 1
                c_sum.TreeOLH += arraytree[l]->children[i]->level_olh;
#endif
            }

            for (int i = 0; i < b; i++)
            {
                if (!arraytree[l]->children[i])
                    break;
                arraytree[l]->children[i]->pert_freq_irr += ((arraytree[l]->pert_freq_irr - c_sum.TreeRR) / c_sum.length);
                arraytree[l]->children[i]->level_ht += ((arraytree[l]->level_ht - c_sum.TreeHT) / c_sum.length);
#if OLH == 1
                arraytree[l]->children[i]->level_olh += ((arraytree[l]->level_olh - c_sum.TreeOLH) / c_sum.length);
#endif
            }

            l++;
        }
    }
}

/*
These is an iterative version of the weighted averaging step.
*/
void B_ary_tree_min::weighted_averaging_itr(node_b_l *root)
{
    int l = 0;
    int r = 0;
    double root_weight = 0.0;
    double children_weight = 0.0;
    for (level = height - 1; level >= 0; level--)
    {
        j = 0;
        l = 0;
        while (j < level)
            l += pow(b, j++);
        r = l + pow(b, level);
        while (l < r)
        {
            root_weight = (pow(b, height - level) - pow(b, height - level - 1.0)) / (pow(b, height - level) - 1.0);
            children_weight = 1.0 - root_weight;
            c_sum.clear();
            for (int i = 0; i < b; i++)
            {
                if (!arraytree[l]->children[i])
                    break;
                c_sum.TreeRR += arraytree[l]->children[i]->pert_freq_irr;
                c_sum.TreeHT += arraytree[l]->children[i]->level_ht;
#if OLH == 1
                c_sum.TreeOLH += arraytree[l]->children[i]->level_olh;
#endif
            }
            arraytree[l]->pert_freq_irr = (arraytree[l]->pert_freq_irr * root_weight) + (c_sum.TreeRR * children_weight);
            arraytree[l]->level_ht = (arraytree[l]->level_ht * root_weight) + (c_sum.TreeHT * children_weight);
#if OLH == 1
            arraytree[l]->level_olh = (arraytree[l]->level_olh * root_weight) + (c_sum.TreeOLH * children_weight);
#endif

            l++;
        }
    }
}

/*
These is a recursive version of the weighted averaging step.
*/
void B_ary_tree_min::weighted_averaging(node_b_l *root)
{
    if (!root)
        return;
    for (int i = 0; i < b; i++)
        weighted_averaging(root->children[i]);

    if ((height - root->level) > 0)
    {
        double root_weight = (pow(b, height - root->level) - pow(b, height - root->level - 1.0)) / (pow(b, height - root->level) - 1.0);
        double children_weight = 1.0 - root_weight;
        //std::cout  << height <<" , "<< root->level << " , " << root->range;
        c_sum.clear();
        for (int i = 0; i < b; i++)
        {
            if (root->children[i])
            {

                c_sum.TreeRR += root->children[i]->pert_freq_irr;
                c_sum.TreeHT += root->children[i]->level_ht;
#if OLH == 1
                c_sum.TreeOLH += root->children[i]->level_olh;
#endif
            }
        }

        root->pert_freq_irr = (root->pert_freq_irr * root_weight) + (c_sum.TreeRR * children_weight);
        root->level_ht = (root->level_ht * root_weight) + (c_sum.TreeHT * children_weight);
#if olh == 1
        root->level_olh = (root->level_olh * root_weight) + (c_sum.TreeOLH * children_weight);
#endif
    }
}

/*
These is a recursive version of the mean consistenct step.
*/

void B_ary_tree_min::mean_Consistency(node_b_l *root)
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
#if OLH == 1
            c_sum.TreeOLH += root->children[i]->level_olh;
#endif
        }
    }

    for (int i = 0; i < b; i++)
    {
        if (root->children[i])
        {
            //cout << abs(-c_sum.sm _rr_u + root->pert_freq_irr_u) / c_sum.length << endl;

            root->children[i]->pert_freq_irr += ((root->pert_freq_irr - c_sum.TreeRR) / c_sum.length);
            root->children[i]->level_ht += ((root->level_ht - c_sum.TreeHT) / c_sum.length);
#if OLH == 1
            root->children[i]->level_olh += ((root->level_olh - c_sum.TreeOLH) / c_sum.length);
#endif
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