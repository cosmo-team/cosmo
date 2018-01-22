#ifndef PACK_COLOR_HPP
#define PACK_COLOR_HPP

int getMilliCount();
int getMilliSpan(int nTimeStart);
typedef struct p
{
    std::string plan_filename = "";
    std::string matrix1_filename = "";
    int num_colors1;
    std::string matrix2_filename = "";
    int num_colors2;
    std::string output_prefix = "";
} parameters_t;


void parse_arguments(int argc, char **argv, parameters_t & params);

// class SDIter {
//     /* from sd_vector, we require:
//        low.size() 
//        high.operator[]
//        low.operator[]
//        wl

//      */
//     sdsl::sd_vector<> *m_s;

//     // position in m_s->low
//     size_t lowndx = 0;

//     // position in m_s->high
//     // we compute a 'next' value during peek() which actually gets locked in by advance()    
//     size_t highndx = 0; // points to next chunk of bits to read from m_s->high
//     size_t nexthigh = 0; // if peek() has been called, points to the next chunk of bits to read from m_s->high

//     // each chunk of bits in m_s->high stores a /delta/.  We need to keep track of the accumulated value of all the deltas
//     // again, we compute a 'next' value during peek() which actually gets locked in by advance()
//     size_t cur_high = 0; // accumulated (i.e. non-delta) value of the position of the last 1 instance
//     size_t next_cur_high = 0; // if peek() has been called, accumulated (i.e. non-delta) value of the position of the peeked 1 instance
//     bool peeked = false;
// public:
//     SDIter(sdsl::sd_vector<> *a_s) {m_s = a_s;};
//     signed long long peek() {
//         if (lowndx >= m_s->low.size()) {
//             return -1;
//         }
//         nexthigh = highndx;
//         next_cur_high = cur_high;
        
//         size_t high_increment = 0; // accumulates number of 0's
//         while (m_s->high[nexthigh] != 1) {
//             high_increment++;
//             nexthigh++;
//         }
//         nexthigh++; // advance to next one
//         next_cur_high += high_increment;
//         size_t low = m_s->low[lowndx];

//         peeked = true;
//         size_t one_loc = (next_cur_high << m_s->wl) | low; // '1' location
//         assert((*m_s)[one_loc]);
//         return one_loc;
//     };
//     void advance() {
//             if (!peeked) { peek();}
//             lowndx++;
//             highndx = nexthigh;
//             cur_high = next_cur_high;
//             peeked = false;

//     }
// };

#endif
