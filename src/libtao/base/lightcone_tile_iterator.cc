#include "lightcone_tile_iterator.hh"

namespace tao
{

lightcone_tile_iterator::lightcone_tile_iterator()
    : _lc(0)
    , _done(true)
    , _idx(std::numeric_limits<unsigned>::max())
{
}

lightcone_tile_iterator::lightcone_tile_iterator(tao::lightcone const& lc)
    : _lc(&lc)
    , _done(false)
    , _idx(0)
{
    LOGBLOCKD("Beginning lightcone tile iterator.");

    // Locate the first box to use. Start by finding a close corner
    // of the lightcone.
    real_type min_ra = lc.min_ra() + lc.ra_offset();
    real_type min_dec = lc.min_dec() + lc.dec_offset();

    std::array<real_type, 3> crd;
    hpc::num::ecs_to_cartesian<real_type>(min_ra, min_dec, crd[0], crd[1], crd[2], lc.min_dist());
    for (int ii = 0; ii < 3; ++ii)
        crd[ii] += lc.origin()[ii];

    LOGDLN("First lightcone corner: ", crd);

    // Make the box corner line up.
    std::array<int, 3> first;
    real_type bs = _lc->simulation()->box_size();
    for (unsigned ii = 0; ii < 3; ++ii)
    {
        if (crd[ii] >= 0.0)
        {
            first[ii] = crd[ii] / bs;
            crd[ii] = first[ii] * bs;
        }
        else
        {
            first[ii] = crd[ii] / bs - 1;
            crd[ii] = first[ii] * bs;
        }
    }

    LOGDLN("Line up corner: ", crd);
    LOGDLN("First index guess: ", first);

    // Check all adjacent boxes, one will be in overlap.
    bool ok = lc.overlap(crd, std::array<real_type, 3>{crd[0] + bs, crd[1] + bs, crd[2] + bs});
    if (!ok)
    {
        first[0]--;
        crd[0] -= bs;
        ok = lc.overlap(crd, std::array<real_type, 3>{crd[0] + bs, crd[1] + bs, crd[2] + bs});
    }
    if (!ok)
    {
        first[1]--;
        crd[1] -= bs;
        ok = lc.overlap(crd, std::array<real_type, 3>{crd[0] + bs, crd[1] + bs, crd[2] + bs});
    }
    if (!ok)
    {
        first[0]++;
        crd[0] += bs;
        ok = lc.overlap(crd, std::array<real_type, 3>{crd[0] + bs, crd[1] + bs, crd[2] + bs});
    }
    if (!ok)
    {
        first[2]--;
        crd[2] -= bs;
        ok = lc.overlap(crd, std::array<real_type, 3>{crd[0] + bs, crd[1] + bs, crd[2] + bs});
    }
    if (!ok)
    {
        first[1]++;
        crd[1] += bs;
        ok = lc.overlap(crd, std::array<real_type, 3>{crd[0] + bs, crd[1] + bs, crd[2] + bs});
    }
    if (!ok)
    {
        first[0]--;
        crd[0] -= bs;
        ok = lc.overlap(crd, std::array<real_type, 3>{crd[0] + bs, crd[1] + bs, crd[2] + bs});
    }
    if (!ok)
    {
        first[1]--;
        crd[1] -= bs;
        ok = lc.overlap(crd, std::array<real_type, 3>{crd[0] + bs, crd[1] + bs, crd[2] + bs});
    }
    ASSERT(ok, "Failed to locate initial overlap.");

    LOGDLN("First corner: ", crd);
    LOGDLN("First index: ", first);

    // Push it onto the stack to be processed, then set it
    // up to be returned.
    _rem_tiles.push_back(first);
    increment();
    _idx = 0;
}

tao::lightcone const* lightcone_tile_iterator::lightcone() const { return _lc; }

bool lightcone_tile_iterator::done() const { return _done; }

unsigned lightcone_tile_iterator::index() const { return _idx; }

std::list<std::array<int, 3>> const& lightcone_tile_iterator::remaining_tiles() const
{
    return _rem_tiles;
}

std::map<std::array<int, 3>, int> const& lightcone_tile_iterator::done_tiles() const
{
    return _done_tiles;
}

void lightcone_tile_iterator::save_checkpoint(boost::property_tree::ptree& pt) const
{
    pt.put("tile", std::to_string(this->index()));
}

void lightcone_tile_iterator::increment()
{
    LOGBLOCKD("Incrementing lightcone tile iterator.");

    // Grab the next tile or finish up.
    if (!_rem_tiles.empty())
    {
        // Increment the index here. I do this because it's a little
        // easier than putting the increment on every return.
        ++_idx;

        // Set the current tile.
        auto pos = _rem_tiles.front();
        _rem_tiles.pop_front();
        real_type bs = _lc->simulation()->box_size();
        std::array<real_type, 3> crd{pos[0] * bs, pos[1] * bs, pos[2] * bs};
        _tile = tile<real_type>(_lc, crd);
        LOGBLOCKD("Setting current tile to: ", crd);

        // Add the current tile to the done set.
        auto cur_done_it = _done_tiles.insert(std::make_pair(pos, 0)).first;
        LOGDLN("Added to done set, which has value: ", cur_done_it->second);

        // Add neighboring tiles to queue.
        for (unsigned ii = 0; ii < 3; ++ii)
        {
            for (int jj = -1; jj < 2; jj += 2)
            {
                pos[ii] += jj;
                LOGBLOCKD("Checking neighbor position: ", pos);
                auto it = _done_tiles.find(pos);
                if (it != _done_tiles.end())
                {
                    LOGDLN("Already in done set with value: ", it->second + 1);
                    if (++it->second == 6)
                    {
                        _done_tiles.erase(it);
                        LOGDLN("Erased.");
                    }
                }
                else
                {
                    LOGDLN("Not in done set.");
                    std::array<real_type, 3> min, max;
                    for (int kk = 0; kk < 3; ++kk)
                    {
                        min[kk] = pos[kk] * bs;
                        max[kk] = min[kk] + bs;
                    }
                    if (_lc->overlap(min, max))
                    {
                        LOGDLN("In overlap, adding to remaining tiles.");
                        _rem_tiles.push_back(pos);
                        _done_tiles.insert(std::make_pair(pos, 1));
                    }
                    else
                    {
                        LOGDLN("Not in overlap, incrementing current done iterator to: ",
                               cur_done_it->second + 1);
                        if (++(cur_done_it->second) == 6)
                        {
                            _done_tiles.erase(cur_done_it);
                            LOGDLN("Erased.");
                        }
                    }
                }
                pos[ii] -= jj;
            }
        }

        LOGDLN("Active tiles: ", _rem_tiles);
        LOGDLN("Done tiles: ", _done_tiles);
    }
    else
    {
        _done = true;
        LOGDLN("Lightcone tile iterator done.");
    }
}

bool lightcone_tile_iterator::equal(const lightcone_tile_iterator& op) const
{
    return _done == op._done;
}

lightcone_tile_iterator::reference_type lightcone_tile_iterator::dereference() const
{
    return _tile;
}

} // namespace tao
