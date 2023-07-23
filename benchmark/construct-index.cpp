#include <iostream>
#include <algorithm>
#include <cmath>

#include <sdsl/config.hpp>
#include <sdsl/construct.hpp>

#include <ri/bwt.h>
#include <ri/rle_string.hpp>
#include "definitions.h"

void construct_text(std::string input_file, sdsl::cache_config *t_config)
{

	std::size_t n;
	sdsl::int_vector<8> text;

	std::string input;
	{
	  std::ifstream fs(input_file);
	  std::stringstream buffer;
	  buffer << fs.rdbuf();
	  
	  input = buffer.str();
	}

	// Construct text representation for SDSL use.
	n = input.size();
	text.resize(input.size() + 1);

	replace_copy(input.begin(), input.end(), text.begin(), 0, 2);

	text[text.size() - 1] = 0; // Append symbol zero at the end

	sdsl::store_to_cache(text, sdsl::conf::KEY_TEXT, *t_config);

}

void construct_sa(sdsl::cache_config *t_config)
{
	// Use SDSL functionality to build the SA
    sdsl::construct_sa<8>(*t_config);
}

void construct_bwt(sdsl::cache_config *t_config)
{
	typedef sdsl::int_vector<>::size_type size_type;
	typedef sdsl::int_vector<8> text_type;
	typedef sdsl::int_vector_buffer<8> bwt_type;

	// Load text from disk
	sdsl::int_vector<8> text;
	load_from_cache(text, sdsl::conf::KEY_TEXT, *t_config);
	auto n = text.size();
	uint8_t bwt_width = text.width();

	sdsl::int_vector<> bwt_heads_pos; // BWT run heads positions in BWT Array
	sdsl::int_vector<> bwt_heads_text_pos; // BWT run heads positions in Text
	std::vector<std::size_t> bwt_heads_text_pos_vec; // BWT run heads positions in Text (vector)
	sdsl::int_vector<> bwt_tails_pos; // BWT run tails positions in BWT Array
	sdsl::int_vector<> bwt_tails_text_pos; // BWT run tails positions in Text
	std::vector<std::size_t> bwt_tails_text_pos_vec; // BWT run tails positions in Text (vector)
	sdsl::bit_vector bwt_heads_in_text_bv; // Mark positions of BWT run heads in text (bitvector)
	sdsl::sd_vector<> bwt_heads_in_text_bv_sd; // Mark positions of BWT run heads in text (sparse bitvector)
	std::vector<std::size_t> f(256, 0); // F Array

	std::size_t r;

	{
		// Prepare to stream SA from disc
		size_type buffer_size = 1000000; // buffer_size is a multiple of 8!, TODO: still true?
		sdsl::int_vector_buffer<> sa_buf(sdsl::cache_file_name(sdsl::conf::KEY_SA, *t_config), std::ios::in, buffer_size);

		// Build BWT sequentially by streaming SA and random access to text
		std::string bwt_file = cache_file_name(sdsl::conf::KEY_BWT, *t_config);
		bwt_type bwt_buf(bwt_file, std::ios::out, buffer_size, bwt_width);

		f.clear();
		f.resize(256, 0);

		auto report_bwt = [&bwt_buf, &f](auto idx, auto symbol) {
		  bwt_buf.push_back(symbol);
		  ++f[symbol + 1];
		};

		// Report BWT run heads
		std::vector<std::size_t> bwt_heads_pos_vec;
		bwt_heads_text_pos_vec.clear();
		auto report_bwt_head =
		    [&bwt_heads_pos_vec, &bwt_heads_text_pos_vec](
		        const auto &run, const auto &idx, const auto &symbol, const auto &pos) {
		      bwt_heads_pos_vec.emplace_back(idx);
		      bwt_heads_text_pos_vec.emplace_back(pos);
		    };

		// Report BWT run heads
		std::vector<std::size_t> bwt_tails_pos_vec;
		bwt_tails_text_pos_vec.clear();
		auto report_bwt_tail =
		    [&bwt_tails_pos_vec, &bwt_tails_text_pos_vec](
		        const auto &run, const auto &idx, const auto &symbol, const auto &pos) {
		      bwt_tails_pos_vec.emplace_back(idx);
		      bwt_tails_text_pos_vec.emplace_back(pos);
		    };

		// Compute BWT and its runs
		r = ri::computeBWT(text.size(),
		                   [&sa_buf](auto idx) { return sa_buf[idx]; },
		                   [&text](auto idx) { return text[idx]; },
		                   report_bwt,
		                   report_bwt_head,
		                   report_bwt_tail);

		assert(text.size() == bwt_buf.size());
		assert(r == bwt_heads_pos_vec.size());
		assert(r == bwt_tails_pos_vec.size());

		auto log_n = sdsl::bits::hi(n) + 1;

		// Build compact representations for heads and tails of the BWT runs
		bwt_heads_pos = sdsl::int_vector<>(r, 0, log_n);
		bwt_heads_text_pos = sdsl::int_vector<>(r, 0, log_n);

		bwt_tails_pos = sdsl::int_vector<>(r, 0, log_n);
		bwt_tails_text_pos = sdsl::int_vector<>(r, 0, log_n);

		bwt_heads_in_text_bv = sdsl::bit_vector(text.size(), 0);

		for (std::size_t i = 0; i < r; ++i) {
		  bwt_heads_pos[i] = bwt_heads_pos_vec[i];
		  bwt_heads_text_pos[i] = bwt_heads_text_pos_vec[i];
		  bwt_heads_in_text_bv[bwt_heads_text_pos_vec[i]] = 1;

		  bwt_tails_pos[i] = bwt_tails_pos_vec[i];
		  bwt_tails_text_pos[i] = bwt_tails_text_pos_vec[i];
		}
		bwt_heads_in_text_bv_sd = sdsl::sd_vector<>(bwt_heads_in_text_bv);

		// Build F array.
		for (std::size_t i = 1; i < 256; ++i) {
		  f[i] += f[i - 1];
		}

		assert(f[255] == n);
	}

  sdsl::store_to_cache(bwt_heads_pos, ri::KEY_BWT_HEADS, *t_config);
  sdsl::store_to_cache(bwt_heads_text_pos, ri::KEY_BWT_HEADS_TEXT_POS, *t_config);
  sdsl::store_to_cache(bwt_heads_text_pos_vec, ri::KEY_BWT_HEADS_TEXT_POS + "_vec", *t_config);
  sdsl::store_to_cache(bwt_heads_in_text_bv, ri::KEY_BWT_HEADS_TEXT_POS + "_bv", *t_config);
  sdsl::store_to_cache(bwt_heads_in_text_bv_sd, ri::KEY_BWT_HEADS_TEXT_POS + "_bv_sd", *t_config);

  sdsl::store_to_cache(bwt_tails_pos, ri::KEY_BWT_TAILS, *t_config);
  sdsl::store_to_cache(bwt_tails_text_pos, ri::KEY_BWT_TAILS_TEXT_POS, *t_config);
  sdsl::store_to_cache(bwt_tails_text_pos_vec, ri::KEY_BWT_TAILS_TEXT_POS + "_vec", *t_config);

  sdsl::store_to_cache(f, ri::KEY_F, *t_config);

}

/// Build run length codification of Burrows Wheeler Transforms.
void construct_bwt_rle(sdsl::cache_config *t_config){
  std::size_t buffer_size = 1000000; // buffer_size is a multiple of 8!, TODO: still true?
  std::string bwt_file = cache_file_name(sdsl::conf::KEY_BWT, *t_config);
  sdsl::int_vector_buffer<8> bwt_buf(bwt_file, std::ios::in, buffer_size);

  std::string bwt_s;
  replace_copy(bwt_buf.begin(), bwt_buf.end(), back_inserter(bwt_s), 0, 1);

  ri::rle_string<> bwt_rle;
  {
    bwt_rle = ri::rle_string<>(bwt_s);
  }

  sdsl::store_to_cache(bwt_rle, ri::KEY_BWT_RLE, *t_config);
};

/// Sort BWT run tails by its positions in the text.
void sort_bwt_tails(sdsl::cache_config *t_config){
  std::vector<std::size_t> bwt_tails_text_pos_vec;
  sdsl::load_from_cache(bwt_tails_text_pos_vec, ri::KEY_BWT_TAILS_TEXT_POS + "_vec", *t_config);

  std::vector<std::size_t> tails_idxs; // Indices of the run tails sorted by its text positions

  {
    tails_idxs.clear();
    tails_idxs.resize(bwt_tails_text_pos_vec.size());
    iota(tails_idxs.begin(), tails_idxs.end(), 0);

    sort(tails_idxs.begin(),
         tails_idxs.end(),
         [&bwt_tails_text_pos_vec](const auto &a, const auto &b) -> bool {
           return bwt_tails_text_pos_vec[a] < bwt_tails_text_pos_vec[b];
         });
  }

  sdsl::store_to_cache(tails_idxs, ri::KEY_BWT_TAILS_TEXT_POS_SORTED_IDX + "_vec", *t_config);
};

/// Sort BWT run heads by its positions in the text.
void sort_bwt_heads(sdsl::cache_config *t_config){
  std::vector<std::size_t> bwt_heads_text_pos_vec;
  sdsl::load_from_cache(bwt_heads_text_pos_vec, ri::KEY_BWT_HEADS_TEXT_POS + "_vec", *t_config);
  auto r = bwt_heads_text_pos_vec.size();
  auto log_r = sdsl::bits::hi(r) + 1;

  std::vector<std::size_t> heads_idxs; // Indices of the run tails sorted by its text positions
  sdsl::int_vector<> tail_idxs_by_heads_in_text;

  {
    heads_idxs.clear();
    heads_idxs.resize(r);
    iota(heads_idxs.begin(), heads_idxs.end(), 0);

    sort(heads_idxs.begin(),
         heads_idxs.end(),
         [&bwt_heads_text_pos_vec](const auto &a, const auto &b) -> bool {
           return bwt_heads_text_pos_vec[a] < bwt_heads_text_pos_vec[b];
         });

    tail_idxs_by_heads_in_text = sdsl::int_vector<>(r, 0, log_r);
    transform(heads_idxs.begin(),
              heads_idxs.end(),
              tail_idxs_by_heads_in_text.begin(),
              [r](auto i) { return (i + r - 1) % r; });
  }

  sdsl::store_to_cache(heads_idxs, ri::KEY_BWT_HEADS_TEXT_POS_SORTED_IDX + "_vec", *t_config);
  sdsl::store_to_cache(tail_idxs_by_heads_in_text,
                       ri::KEY_BWT_TAILS_SAMPLED_IDX_BY_HEAD_IN_TEXT,
                       *t_config);
};

/// Build BWT run tails sampling.
void tails_sampling(sdsl::cache_config *t_config, std::size_t s){

  std::vector<std::size_t> bwt_tails_text_pos_vec;
  sdsl::load_from_cache(bwt_tails_text_pos_vec, ri::KEY_BWT_TAILS_TEXT_POS + "_vec", *t_config);

  std::vector<std::size_t> tails_idxs_sorted;
  sdsl::load_from_cache(tails_idxs_sorted, ri::KEY_BWT_TAILS_TEXT_POS_SORTED_IDX + "_vec", *t_config);
  std::size_t r = tails_idxs_sorted.size();

  // We must sample the run-tails previous to the first and last run-head in the text
  std::size_t prev_tail_to_first_head_in_text;
  std::size_t prev_tail_to_last_head_in_text;
  {
    std::vector<std::size_t> bwt_head_idxs_sorted_in_text;
    sdsl::load_from_cache(bwt_head_idxs_sorted_in_text, ri::KEY_BWT_HEADS_TEXT_POS_SORTED_IDX + "_vec", *t_config);

    assert(bwt_head_idxs_sorted_in_text.size() == r);

    prev_tail_to_first_head_in_text = (bwt_head_idxs_sorted_in_text.front() + r - 1) % r;
    prev_tail_to_last_head_in_text = (bwt_head_idxs_sorted_in_text.back() + r - 1) % r;
  }

  std::vector<std::size_t> sampled_idxs_vec; // Indices of sampled BWT run tails
  sdsl::int_vector<> sampled_tails_in_text; // Text position of sampled BWT run tails in order of the BWT Array
  sdsl::bit_vector sampled_tails_idx_bv; // Mark which tails are sampled (bitvector)
  sdsl::sd_vector<> sampled_tails_idx_bv_sd; // Mark which tails are sampled (sparse bitvector)

  std::size_t r_prime = 0;

  {
    sampled_idxs_vec.clear();
    sampled_idxs_vec.reserve(tails_idxs_sorted.size() / 2);

    // Compute sampled BWT run tails
    auto last_sampled_idx = tails_idxs_sorted[0];
    auto last_sampled = bwt_tails_text_pos_vec[last_sampled_idx];
    sampled_idxs_vec.emplace_back(last_sampled_idx);

    std::size_t prev_idx = tails_idxs_sorted[1];
    std::size_t prev = bwt_tails_text_pos_vec[prev_idx];
    for (std::size_t i = 2; i < tails_idxs_sorted.size(); ++i) {
      auto current_idx = tails_idxs_sorted[i];
      auto current = bwt_tails_text_pos_vec[current_idx];

      if (prev_idx == prev_tail_to_first_head_in_text ||
          prev_idx == prev_tail_to_last_head_in_text ||
          s < current - last_sampled) {
        // Sample the previous run tails in text
        sampled_idxs_vec.emplace_back(prev_idx);
        last_sampled = prev;
      }

      prev_idx = current_idx;
      prev = current;
    }

    sampled_idxs_vec.emplace_back(prev_idx);
    r_prime = sampled_idxs_vec.size();

    // Sort indexes of sampled BWT run tails in order of the BWT Array
    sort(sampled_idxs_vec.begin(), sampled_idxs_vec.end());

    auto log_n = sdsl::bits::hi(prev) + 1;
    sampled_tails_in_text = sdsl::int_vector<>(r_prime, 0, log_n);

    sampled_tails_idx_bv = sdsl::bit_vector(r, 0);

    // Compute text position of the sampled BWT run tails and mark the sampled tails indices in a bitvector
    for (std::size_t i = 0; i < r_prime; ++i) {
      sampled_tails_in_text[i] = bwt_tails_text_pos_vec[sampled_idxs_vec[i]];
      sampled_tails_idx_bv[sampled_idxs_vec[i]] = 1;
    }

    sampled_tails_idx_bv_sd = sdsl::sd_vector<>(sampled_tails_idx_bv);
  }

  sdsl::store_to_cache(sampled_tails_in_text, std::to_string(s) + "_" + ri::KEY_BWT_TAILS_TEXT_POS_SAMPLED, *t_config);
  sdsl::store_to_cache(sampled_idxs_vec, std::to_string(s) + "_" + ri::KEY_BWT_TAILS_SAMPLED_IDX + "_vec", *t_config);
  sdsl::store_to_cache(sampled_tails_idx_bv,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_SAMPLED_IDX + "_bv",
                       *t_config);
  sdsl::store_to_cache(sampled_tails_idx_bv_sd,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_SAMPLED_IDX + "_bv_sd",
                       *t_config);
};

/// Build BWT run heads sampling.
void heads_sampling(sdsl::cache_config *t_config, std::size_t s){

  std::vector<std::size_t> heads_in_text_vec;
  sdsl::load_from_cache(heads_in_text_vec, ri::KEY_BWT_HEADS_TEXT_POS + "_vec", *t_config);
  std::size_t r = heads_in_text_vec.size();

  std::vector<std::size_t> sampled_tail_idxs_vec;
  sdsl::load_from_cache(sampled_tail_idxs_vec,
                        std::to_string(s) + "_" + ri::KEY_BWT_TAILS_SAMPLED_IDX + "_vec",
                        *t_config);
  std::size_t r_prime = sampled_tail_idxs_vec.size();
  auto log_r_prime = sdsl::bits::hi(r_prime) + 1;

  std::vector<std::size_t> head_idxs_sorted_in_text;
  sdsl::load_from_cache(head_idxs_sorted_in_text, ri::KEY_BWT_HEADS_TEXT_POS_SORTED_IDX + "_vec", *t_config);
  assert(head_idxs_sorted_in_text.size() == r);

  // Indices of sampled BWT run tails indices, sorted by the text position of its following run heads (vector)
  std::vector<std::size_t> sampled_tail_idx_by_heads_in_text_vec;
  // Indices of sampled BWT run tails indices, sorted by the text position of its following run heads
  sdsl::int_vector<> sampled_tail_idx_by_heads_in_text;
  // Marked trustworthy indices of the sampled BWT run heads: bv[i] = true iff between sampled run heads i and i + 1 none run head was deleted
  std::vector<bool> marked_sampled_idxs_vec;
  sdsl::bit_vector marked_sampled_idxs_bv;
  sdsl::sd_vector<> marked_sampled_idxs_bv_sd;
  // Trusted interval for marked trustworthy indices of the sampled BWT run heads
  std::vector<std::size_t> trusted_area_for_marked_sampled_idxs_vec;
  sdsl::int_vector<> trusted_area_for_marked_sampled_idxs;
  sdsl::bit_vector sampled_heads_in_text_bv; // Mark positions of sampled BWT run heads in text (bitvector)
  sdsl::sd_vector<> sampled_heads_in_text_bv_sd; // Mark positions of sampled BWT run heads in text (sparse bitvector)

  {
    sampled_tail_idx_by_heads_in_text_vec.clear();
    sampled_tail_idx_by_heads_in_text_vec.resize(r_prime);
    std::iota(sampled_tail_idx_by_heads_in_text_vec.begin(), sampled_tail_idx_by_heads_in_text_vec.end(), 0);

    // Sort indices of sampled BWT run tails indices by text position of the next run heads in BWT Array
    sort(sampled_tail_idx_by_heads_in_text_vec.begin(),
         sampled_tail_idx_by_heads_in_text_vec.end(),
         [&heads_in_text_vec, &sampled_tail_idxs_vec, r](const auto &a, const auto &b) {
           auto idx_a = (sampled_tail_idxs_vec[a] + 1) % r;
           auto idx_b = (sampled_tail_idxs_vec[b] + 1) % r;
           return heads_in_text_vec[idx_a] < heads_in_text_vec[idx_b];
         });

    assert((sampled_tail_idxs_vec[sampled_tail_idx_by_heads_in_text_vec.front()] + 1) % r
               == head_idxs_sorted_in_text.front());
    assert((sampled_tail_idxs_vec[sampled_tail_idx_by_heads_in_text_vec.back()] + 1) % r
               == head_idxs_sorted_in_text.back());

    marked_sampled_idxs_vec.clear();
    marked_sampled_idxs_vec.resize(r_prime, true);
    marked_sampled_idxs_bv = sdsl::bit_vector(r_prime, true);

    auto n = heads_in_text_vec[head_idxs_sorted_in_text.back()] + 1;
    sampled_heads_in_text_bv = sdsl::bit_vector(n, 0);

    trusted_area_for_marked_sampled_idxs_vec.clear();
    trusted_area_for_marked_sampled_idxs_vec.reserve(r_prime / 2);

    // Marks sampled BWT run heads indices in text and if they are trustworthy
    {
      auto sampled_head_idx = head_idxs_sorted_in_text.front();
      for (std::size_t i = 0, j = 0; i < r_prime - 1; ++i) {
        auto sampled_head_tpos = heads_in_text_vec[sampled_head_idx];
        // Marks sampled BWT run heads indices in text
        sampled_heads_in_text_bv[sampled_head_tpos] = true;

        auto next_sampled_head_idx = (sampled_tail_idxs_vec[sampled_tail_idx_by_heads_in_text_vec[i + 1]] + 1) % r;
        auto head_idx = head_idxs_sorted_in_text[++j];
        if (next_sampled_head_idx != head_idx) {
          // Unmarks untrusted sampled BWT run heads indices
          marked_sampled_idxs_vec[i] = false;
          marked_sampled_idxs_bv[i] = false;

          // Store trusted range for untrusted sampled BWT run heads indices
          trusted_area_for_marked_sampled_idxs_vec.emplace_back(heads_in_text_vec[head_idx] - sampled_head_tpos + 1);

          while (next_sampled_head_idx != head_idxs_sorted_in_text[++j]);
        }

        sampled_head_idx = next_sampled_head_idx;
      }

      sampled_heads_in_text_bv[heads_in_text_vec[head_idxs_sorted_in_text.back()]] = true;
    }
    sampled_heads_in_text_bv_sd = sdsl::sd_vector<>(sampled_heads_in_text_bv);
    marked_sampled_idxs_bv_sd = sdsl::sd_vector<>(marked_sampled_idxs_bv);

    sampled_tail_idx_by_heads_in_text = sdsl::int_vector<>(r_prime, 0, log_r_prime);
    copy(sampled_tail_idx_by_heads_in_text_vec.begin(),
         sampled_tail_idx_by_heads_in_text_vec.end(),
         sampled_tail_idx_by_heads_in_text.begin());

    trusted_area_for_marked_sampled_idxs = sdsl::int_vector<>(trusted_area_for_marked_sampled_idxs_vec.size(), 0);
    copy(trusted_area_for_marked_sampled_idxs_vec.begin(),
         trusted_area_for_marked_sampled_idxs_vec.end(),
         trusted_area_for_marked_sampled_idxs.begin());
    sdsl::util::bit_compress(trusted_area_for_marked_sampled_idxs);
  }

  sdsl::store_to_cache(sampled_tail_idx_by_heads_in_text,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_SAMPLED_IDX_BY_HEAD_IN_TEXT,
                       *t_config);
  sdsl::store_to_cache(marked_sampled_idxs_vec,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_MARKED_SAMPLED_IDX_BY_HEAD_IN_TEXT + "_vec",
                       *t_config);
  sdsl::store_to_cache(marked_sampled_idxs_bv,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_MARKED_SAMPLED_IDX_BY_HEAD_IN_TEXT + "_bv",
                       *t_config);
  sdsl::store_to_cache(marked_sampled_idxs_bv_sd,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_MARKED_SAMPLED_IDX_BY_HEAD_IN_TEXT + "_bv_sd",
                       *t_config);
  sdsl::store_to_cache(sampled_heads_in_text_bv,
                       std::to_string(s) + "_" + ri::KEY_BWT_HEADS_SAMPLED_TEXT_POS + "_bv",
                       *t_config);
  sdsl::store_to_cache(sampled_heads_in_text_bv_sd,
                       std::to_string(s) + "_" + ri::KEY_BWT_HEADS_SAMPLED_TEXT_POS + "_bv_sd",
                       *t_config);
  sdsl::store_to_cache(trusted_area_for_marked_sampled_idxs,
                       std::to_string(s) + "_" + ri::KEY_BWT_HEADS_MARKED_SAMPLED_TRUSTED_AREA_IN_TEXT,
                       *t_config);
};



int main(int argc, char **argv) {

		std::string input_file = std::string(argv[1]);
		sdsl::cache_config config(false, ".", sdsl::util::basename(input_file));
    std::cout << "Computing the sr-index of: " << input_file << "\n##########################\n";

    std::cout << "constructing sdsl text representation\n";
		construct_text(input_file, &config);
    std::cout << "constructing the suffix array\n";
		construct_sa(&config);
    std::cout << "constructing the bwt\n";
		construct_bwt(&config);
    std::cout << "constructing the rle bwt\n";
		construct_bwt_rle(&config);
    std::cout << "sorting the bwt tails\n";
    sort_bwt_tails(&config);
    std::cout << "sorting the bwt heads\n";
		sort_bwt_heads(&config);

		for(std::size_t i=0; i<5; ++i)
		{
				std::size_t s = 4 * std::pow(2,i);

        std::cout << "sampling tails with s: " << s << '\n';
				tails_sampling(&config, s);
        std::cout << "sampling heads with s: " << s << '\n';
				heads_sampling(&config, s);
		}

    return 0;
}