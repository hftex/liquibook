// Copyright (c) 2012 - 2017 Object Computing, Inc.
// All rights reserved.
// See the file license.txt for licensing information.
#pragma once

#include "callback.h"
#include "comparable_price.h"
#include "logger.h"
#include "order_book_listener.h"
#include "order_listener.h"
#include "order_tracker.h"
#include "trade_listener.h"
#include "version.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <list>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef LIQUIBOOK_IGNORES_DEPRECATED_CALLS
#define COMPLAIN_ONCE(message)
#else // LIQUIBOOK_IGNORES_DEPRECATED_CALLS
#define COMPLAIN_ONCE(message)                                                 \
  do {                                                                         \
    static bool once = true;                                                   \
    if (once) {                                                                \
      once = false;                                                            \
      std::cerr << "One-time Warning: " << message << std::endl;               \
    }                                                                          \
  } while (false)
#endif // LIQUIBOOK_IGNORES_DEPRECATED_CALLS

namespace liquibook {
namespace book {

struct AcceptTag {};
struct PreAcceptTag {};
struct FillTag {};
struct RejectTag {};
struct CancelTag {};
struct CancelRejectTag {};
struct ReplaceTag {};
struct ReplaceRejectTag {};
struct UpdateTag {};

template <typename OrderPtr> class OrderListener;

template <class OrderBook> class OrderBookListener;

/// @brief The limit order book of a security.  Template implementation allows
///        user to supply common or smart pointers, and to provide a different
///        Order class completely (as long as interface is obeyed).
template <typename OrderPtr,
          template <typename...> typename Multimap = std::multimap>
class OrderBook {
public:
  typedef OrderTracker<OrderPtr> Tracker;
  typedef Callback<OrderPtr, Multimap> TypedCallback;
  typedef OrderListener<OrderPtr> TypedOrderListener;
  typedef OrderBook<OrderPtr, Multimap> MyClass;
  typedef TradeListener<MyClass> TypedTradeListener;
  typedef OrderBookListener<MyClass> TypedOrderBookListener;
  typedef std::vector<TypedCallback> Callbacks;
  typedef Multimap<ComparablePrice, Tracker> TrackerMap;
  typedef std::vector<Tracker> TrackerVec;
  // Keep this around briefly for compatibility.
  typedef TrackerMap Bids;
  typedef TrackerMap Asks;

  typedef std::list<typename TrackerMap::iterator> DeferredMatches;

  /// @brief construct
  OrderBook(const std::string &symbol = "unknown");

  /// @brief Set symbol for orders in this book.
  void set_symbol(const std::string &symbol);

  /// @ Get the symbol for orders in this book
  /// @return the symbol.
  const std::string &symbol() const;

  /// @brief set the order listener
  void set_order_listener(TypedOrderListener *listener);

  /// @brief set the trade listener
  void set_trade_listener(TypedTradeListener *listener);

  /// @brief set the order book listener
  void set_order_book_listener(TypedOrderBookListener *listener);

  /// @brief let the application handle reporting errors.
  void set_logger(Logger *logger);

  /// @brief add an order to book
  /// @param order the order to add
  /// @param conditions special conditions on the order
  /// @return true if the add resulted in a fill
  template <typename Listener>
  bool add(const OrderPtr &order, OrderConditions conditions,
           Listener &&listener);

  /// @brief cancel an order in the book
  template <typename Listener>
  void cancel(const OrderPtr &order, Listener &&listener);

  /// @brief replace an order in the book
  /// @param order the order to replace
  /// @param size_delta the change in size for the order (positive or negative)
  /// @param new_price the new order price, or PRICE_UNCHANGED
  /// @return true if the replace resulted in a fill
  template <typename Listener>
  bool replace(const OrderPtr &order, int32_t size_delta, Price new_price,
               Listener &&listener);

  /// @brief Set the current market price
  /// Intended to be used during initialization to establish the market
  /// price before this order book has generated any exceptions.
  ///
  /// If price is zero (or never set) no market-to-market trades can happen.
  /// @param price is the current market price for this security.
  void set_market_price(Price price);

  /// @brief Get current market price.
  /// The market price is normally the price at which the last trade happened.
  Price market_price() const;

  /// @brief access the bids container
  const TrackerMap &bids() const { return bids_; };

  /// @brief access the asks container
  const TrackerMap &asks() const { return asks_; };

  /// @brief access stop bid orders
  const TrackerMap &stopBids() const { return stopBids_; }

  /// @brief access stop ask orders
  const TrackerMap &stopAsks() const { return stopAsks_; }

  /// @brief move callbacks to another thread's container
  /// @deprecated  This doesn't do anything now
  /// so don't bother to call it in new code.
  void move_callbacks(Callbacks &target);

  /// @brief log the orders in the book.
  std::ostream &log(std::ostream &out) const;

protected:
  /// @brief Internal method to process callbacks.
  /// Protected against recursive calls in case callbacks
  /// issue new requests.
  // void callback_now();

  /// @brief match a new order to current orders
  /// @param inbound_order the inbound order
  /// @param inbound_price price of the inbound order
  /// @param current_orders open orders
  /// @param[OUT] deferred_aons AON orders from current_orders
  ///             that matched the inbound price,
  ///             but were not filled due to quantity
  /// @return true if a match occurred
  template <typename Listener>
  bool match_order(Tracker &inbound_order, Price inbound_price,
                   TrackerMap &current_orders, DeferredMatches &deferred_aons,
                   Listener &&listener);

  template <typename Listener>
  bool match_aon_order(Tracker &inbound, Price inbound_price,
                       TrackerMap &current_orders,
                       DeferredMatches &deferred_aons, Listener &&listener);

  template <typename Listener>
  bool match_regular_order(Tracker &inbound, Price inbound_price,
                           TrackerMap &current_orders,
                           DeferredMatches &deferred_aons, Listener &&listener);

  template <typename Listener>
  Quantity
  try_create_deferred_trades(Tracker &inbound,
                             DeferredMatches &deferred_matches,
                             Quantity maxQty, // do not exceed
                             Quantity minQty, // must be at least
                             TrackerMap &current_orders Listener &&listener);

  /// @brief see if any deferred All Or None orders can now execute.
  /// @param aons iterators to the orders that might now match
  /// @param deferredTrackers the container of the aons
  /// @param marketTrackers the orders to check for matches
  bool check_deferred_aons(DeferredMatches &aons, TrackerMap &deferredTrackers,
                           TrackerMap &marketTrackers);

  /// @brief perform fill on two orders
  /// @param inbound_tracker the new (or changed) order tracker
  /// @param current_tracker the current order tracker
  /// @param max_quantity maximum quantity to trade.
  /// @return the number of units traded (zero if unsuccessful).
  template <typename Listener>
  Quantity create_trade(Tracker &inbound_tracker, Tracker &current_tracker,
                        Quantity max_quantity, Listener &&listener);

  /// @brief find an order in a container
  /// @param order is the the order we are looking for
  /// @param sideMap contains the container where we will look
  /// @param[OUT] result will point to the entry in the container if we find a
  /// match
  /// @returns true: match, false: no match
  bool find_on_market(const OrderPtr &order,
                      typename TrackerMap::iterator &result);

  /// @brief add incoming stop order to stops colletion unless it's already
  /// on the market.
  /// @return true if added to stops, false if it should go directly to the
  /// order book.
  bool add_stop_order(Tracker &tracker);

  /// @brief See if any stop orders should go on the market.
  void check_stop_orders(bool side, Price price, TrackerMap &stops);

  /// @brief accept pending (formerly stop) orders.
  void submit_pending_orders();

  ///////////////////////////////
  // Callback interfaces as
  // virtual methods to simplify
  // derived classes.
  ///////////////////////////////
  // Order Listener interface
  /// @brief callback for an order accept
  virtual void on_accept(const OrderPtr &order, Quantity quantity) {}

  /// @brief callback for an order reject
  virtual void on_reject(const OrderPtr &order, const char *reason) {}

  /// @brief callback for an order fill
  /// @param order the inbound order
  /// @param matched_order the matched order
  /// @param fill_qty the quantity of this fill
  /// @param fill_cost the cost of this fill (qty * price)
  virtual void on_fill(const OrderPtr &order, const OrderPtr &matched_order,
                       Quantity fill_qty, Cost fill_cost,
                       bool inbound_order_filled, bool matched_order_filled) {}

  /// @brief callback for an order cancellation
  virtual void on_cancel(const OrderPtr &order, Quantity quantity) {}

  /// @brief callback for an order cancel rejection
  virtual void on_cancel_reject(const OrderPtr &order, const char *reason) {}

  /// @brief callback for an order replace
  /// @param order the replaced order
  /// @param size_delta the change to order quantity
  /// @param new_price the updated order price
  virtual void on_replace(const OrderPtr &order, Quantity current_qty,
                          Quantity new_qty, Price new_price) {}

  /// @brief callback for an order replace rejection
  virtual void on_replace_reject(const OrderPtr &order, const char *reason) {}

  // End of OrderListener Interface
  ///////////////////////////////
  // TradeListener Interface
  /// @brief callback for a trade
  /// @param book the order book of the fill (not defined whether this is before
  ///      or after fill)
  /// @param qty the quantity of this fill
  /// @param cost the cost of this fill (qty * price)
  virtual void on_trade(const OrderBook *book, Quantity qty, Cost cost) {}
  // End of TradeListener Interface
  ///////////////////////////////
  // BookListener Interface
  /// @brief callback for change anywhere in order book
  virtual void on_order_book_change() {}
  // End of BookListener Interface
  ///////////////////////////////

private:
  bool submit_order(Tracker &inbound);
  bool add_order(Tracker &order_tracker, Price order_price);

private:
  std::string symbol_;
  TrackerMap bids_;
  TrackerMap asks_;

  TrackerMap stopBids_;
  TrackerMap stopAsks_;
  TrackerVec pendingOrders_;

  // Callbacks callbacks_;
  // bool handling_callbacks_;
  TypedOrderListener *order_listener_;
  TypedTradeListener *trade_listener_;
  TypedOrderBookListener *order_book_listener_;
  Logger *logger_;
  Price marketPrice_;
};

template <typename OrderPtr, template <typename...> typename Multimap>
OrderBook<OrderPtr, Multimap>::OrderBook(const std::string &symbol)
    : symbol_(symbol),
      // handling_callbacks_(false),
      order_listener_(nullptr), trade_listener_(nullptr),
      order_book_listener_(nullptr), logger_(nullptr),
      marketPrice_(MARKET_ORDER_PRICE) {
  // callbacks_.reserve(1024); // Why 16?  Why not?
}

template <typename OrderPtr, template <typename...> typename Multimap>
void OrderBook<OrderPtr, Multimap>::set_logger(Logger *logger) {
  logger_ = logger;
}

template <typename OrderPtr, template <typename...> typename Multimap>
void OrderBook<OrderPtr, Multimap>::set_symbol(const std::string &symbol) {
  symbol_ = symbol;
}

template <typename OrderPtr, template <typename...> typename Multimap>
const std::string &OrderBook<OrderPtr, Multimap>::symbol() const {
  return symbol_;
}

template <typename OrderPtr, template <typename...> typename Multimap>
void OrderBook<OrderPtr, Multimap>::set_market_price(Price price) {
  Price oldMarketPrice = marketPrice_;
  marketPrice_ = price;
  if (price > oldMarketPrice || oldMarketPrice == MARKET_ORDER_PRICE) {
    // price has gone up: check stop bids
    bool buySide = true;
    check_stop_orders(buySide, price, stopBids_);
  } else if (price < oldMarketPrice || oldMarketPrice == MARKET_ORDER_PRICE) {
    // price has gone down: check stop asks
    bool buySide = false;
    check_stop_orders(buySide, price, stopAsks_);
  }
}

/// @brief Get current market price.
/// The market price is normally the price at which the last trade happened.
template <typename OrderPtr, template <typename...> typename Multimap>
Price OrderBook<OrderPtr, Multimap>::market_price() const {
  return marketPrice_;
}

template <typename OrderPtr, template <typename...> typename Multimap>
void OrderBook<OrderPtr, Multimap>::set_order_listener(
    TypedOrderListener *listener) {
  order_listener_ = listener;
}

template <typename OrderPtr, template <typename...> typename Multimap>
void OrderBook<OrderPtr, Multimap>::set_trade_listener(
    TypedTradeListener *listener) {
  trade_listener_ = listener;
}

template <typename OrderPtr, template <typename...> typename Multimap>
void OrderBook<OrderPtr, Multimap>::set_order_book_listener(
    TypedOrderBookListener *listener) {
  order_book_listener_ = listener;
}

template <typename OrderPtr, template <typename...> typename Multimap,
          typename Listener>
bool OrderBook<OrderPtr, Multimap>::add(const OrderPtr &order,
                                        OrderConditions conditions,
                                        Listener &&listener) {
  bool matched = false;

  // If the order is invalid, ignore it
  if (order->order_qty() == 0) {
    listener(RejectTag{}, order, "size must be positive");
  } else {
    // size_t accept_cb_index = callbacks_.size();
    listener(PreAcceptTag{}, order);
    Tracker inbound(order, conditions);
    if (inbound.ptr()->stop_price() != 0 && add_stop_order(inbound)) {
      // The order has been added to stops
    } else {
      matched = submit_order(inbound);
      // Note the filled qty in the accept callback
      if (!inbound.filled()) {
        listener(AcceptTag{}, order, inbound.filled_qty());
      }
      // callbacks_[accept_cb_index].quantity = inbound.filled_qty();

      // Cancel any unfilled IOC order
      if (inbound.immediate_or_cancel() && !inbound.filled()) {
        // NOTE - this may need he actual open qty???
        listener(CancelTag{}, order, 0);
      }
    }
    // If adding this order triggered any stops
    // handle those stops now
    while (!pendingOrders_.empty()) {
      submit_pending_orders();
    }
    listener(UpdateTag{});
  }
  // callback_now();
  return matched;
}

template <typename OrderPtr, template <typename...> typename Multimap,
          typename Listener>
void OrderBook<OrderPtr, Multimap>::cancel(const OrderPtr &order,
                                           Listener &&listener) {
  bool found = false;
  Quantity open_qty;
  // If the cancel is a buy order
  if (order->is_buy()) {
    typename TrackerMap::iterator bid;
    find_on_market(order, bid);
    if (bid != bids_.end()) {
      open_qty = bid->second.open_qty();
      // Remove from container for cancel
      bids_.erase(bid);
      found = true;
    }
    // Else the cancel is a sell order
  } else {
    typename TrackerMap::iterator ask;
    find_on_market(order, ask);
    if (ask != asks_.end()) {
      open_qty = ask->second.open_qty();
      // Remove from container for cancel
      asks_.erase(ask);
      found = true;
    }
  }
  // If the cancel was found, issue callback
  if (found) {
    listener(CancelTag{}, order, open_qty);
    listener(Update{});
  } else {
    listener(CancelRejectTag{}, order, "not found");
  }
  // callback_now();
}

template <typename OrderPtr, template <typename...> typename Multimap,
          typename Listener>
bool OrderBook<OrderPtr, Multimap>::replace(const OrderPtr &order,
                                            int32_t size_delta, Price new_price,
                                            Listener &&listener) {
  bool matched = false;
  bool price_change = new_price && (new_price != order->price());

  Price price = (new_price == PRICE_UNCHANGED) ? order->price() : new_price;

  // If the order to replace is a buy order
  TrackerMap &market = order->is_buy() ? bids_ : asks_;
  typename TrackerMap::iterator pos;
  if (find_on_market(order, pos)) {
    // If this is a valid replace
    const Tracker &tracker = pos->second;
    // If there is not enough open quantity for the size reduction
    if (size_delta < 0 && ((int)tracker.open_qty() < -size_delta)) {
      // get rid of as much as we can
      size_delta = -int(tracker.open_qty());
      if (size_delta == 0) {
        // if there is nothing to get rid of
        // Reject the replace
        listener(ReplaceRejectTag{}, tracker.ptr(), "order is already filled");
        return false;
      }
    }

    // Accept the replace
    listener(ReplaceTag{}, order, pos->second.open_qty(), size_delta, price);
    Quantity new_open_qty = pos->second.open_qty() + size_delta;
    pos->second.change_qty(size_delta); // Update my copy
    // If the size change will close the order
    if (!new_open_qty) {
      // Cancel with NO open qty (should be zero after replace)
      listener(CancelTag{}, order, 0);
      market.erase(pos); // Remove order
    } else {
      // Else rematch the new order - there could be a price change
      // or size change - that could cause all or none match
      auto order = pos->second;
      market.erase(pos);                 // Remove old order order
      matched = add_order(order, price); // Add order
    }
    // If replace any order this order triggered any trades
    // which triggered any stops
    // handle those stops now
    while (!pendingOrders_.empty()) {
      submit_pending_orders();
    }
    listener(UpdateTag{});
  } else {
    // not found
    listener(ReplaceRejectTag{}, order, "not found");
  }
  // callback_now();
  return matched;
}

template <typename OrderPtr, template <typename...> typename Multimap>
bool OrderBook<OrderPtr, Multimap>::add_stop_order(Tracker &tracker) {
  bool isBuy = tracker.ptr()->is_buy();
  ComparablePrice key(isBuy, tracker.ptr()->stop_price());
  // if the market price is a better deal then the stop price, it's not time to
  // panic
  bool isStopped = key < marketPrice_;
  if (isStopped) {
    if (isBuy) {
      stopBids_.emplace(key, std::move(tracker));
    } else {
      stopAsks_.emplace(key, std::move(tracker));
    }
  }
  return isStopped;
}

template <typename OrderPtr, template <typename...> typename Multimap>
void OrderBook<OrderPtr, Multimap>::check_stop_orders(bool side, Price price,
                                                      TrackerMap &stops) {
  ComparablePrice until(side, price);
  auto pos = stops.begin();
  while (pos != stops.end()) {
    auto here = pos++;
    if (until < here->first) {
      break;
    }
    pendingOrders_.push_back(std::move(here->second));
    stops.erase(here);
  }
}

template <typename OrderPtr, template <typename...> typename Multimap>
void OrderBook<OrderPtr, Multimap>::submit_pending_orders() {
  TrackerVec pending;
  pending.swap(pendingOrders_);
  for (auto pos = pending.begin(); pos != pending.end(); ++pos) {
    Tracker &tracker = *pos;
    submit_order(tracker);
  }
}

template <typename OrderPtr, template <typename...> typename Multimap>
bool OrderBook<OrderPtr, Multimap>::submit_order(Tracker &inbound) {
  Price order_price = inbound.ptr()->price();
  return add_order(inbound, order_price);
}

template <typename OrderPtr, template <typename...> typename Multimap>
bool OrderBook<OrderPtr, Multimap>::find_on_market(
    const OrderPtr &order, typename TrackerMap::iterator &result) {
  const ComparablePrice key(order->is_buy(), order->price());
  TrackerMap &sideMap = order->is_buy() ? bids_ : asks_;

  for (result = sideMap.find(key); result != sideMap.end(); ++result) {
    // If this is the correct bid
    if (result->second.ptr() == order) {
      return true;
    } else if (key < result->first) {
      // exit early if result is beyond the matching prices
      result = sideMap.end();
      return false;
    }
  }
  return false;
}

// Try to match order.  Generate trades.
// If not completely filled and not IOC,
// add the order to the order book
template <typename OrderPtr, template <typename...> typename Multimap>
bool OrderBook<OrderPtr, Multimap>::add_order(Tracker &inbound,
                                              Price order_price) {
  DeferredMatches deferred_aons;
  // Try to match with current orders
  bool is_buy = inbound.ptr()->is_buy();
  auto &asks = is_buy ? asks_ : bids_;
  auto &bids = is_buy ? bids_ : asks_;
  bool matched = match_order(inbound, order_price, asks, deferred_aons);
  if (!inbound.filled() && !inbound.immediate_or_cancel()) {
    // Insert into 'bids'
    bids.insert({ComparablePrice(is_buy, order_price), inbound});
    // and see if that satisfies any 'ask' orders
    matched |= check_deferred_aons(deferred_aons, asks, bids);
  }
  return matched;
}

template <typename OrderPtr, template <typename...> typename Multimap>
bool OrderBook<OrderPtr, Multimap>::check_deferred_aons(
    DeferredMatches &aons, TrackerMap &deferredTrackers,
    TrackerMap &marketTrackers) {
  bool result = false;
  DeferredMatches ignoredAons;

  for (auto pos = aons.begin(); pos != aons.end(); ++pos) {
    auto entry = *pos;
    ComparablePrice current_price = entry->first;
    Tracker &tracker = entry->second;
    bool matched = match_order(tracker, current_price.price(), marketTrackers,
                               ignoredAons);
    result |= matched;
    if (tracker.filled()) {
      deferredTrackers.erase(entry);
    }
  }
  return result;
}

///  Try to match order at 'price' against 'current' orders
///  If successful
///    generate trade(s)
///    if any current order is complete, remove from 'current' orders
template <typename OrderPtr, template <typename...> typename Multimap>
bool OrderBook<OrderPtr, Multimap>::match_order(
    Tracker &inbound, Price inbound_price, TrackerMap &current_orders,
    DeferredMatches &deferred_aons) {
  if (inbound.all_or_none()) {
    return match_aon_order(inbound, inbound_price, current_orders,
                           deferred_aons);
  }
  return match_regular_order(inbound, inbound_price, current_orders,
                             deferred_aons);
}

template <typename OrderPtr, template <typename...> typename Multimap,
          typename Listener>
bool OrderBook<OrderPtr, Multimap>::match_regular_order(
    Tracker &inbound, Price inbound_price, TrackerMap &current_orders,
    DeferredMatches &deferred_aons, Listener &&listener) {
  // while incoming ! satisfied
  //   current is reg->trade
  //   current is AON:
  //    incoming satisfies AON ->TRADE
  //    add AON to deferred
  // loop
  bool matched = false;
  Quantity inbound_qty = inbound.open_qty();
  // typename TrackerMap::iterator pos = current_orders.begin();
  while (!current_orders.empty() && !inbound.filled()) {
    auto entry = current_orders.begin();
    const ComparablePrice &current_price = entry->first;
    if (!current_price.matches(inbound_price)) {
      // no more trades against current orders are possible
      break;
    }

    //////////////////////////////////////
    // Current price matches inbound price
    Tracker &current_order = entry->second;
    Quantity current_quantity = current_order.open_qty();

    if (current_order.all_or_none()) {
      // if the inbound order can satisfy the current order's AON condition
      if (current_quantity <= inbound_qty) {
        // current is AON, inbound is not AON.
        // inbound can satisfy current's AON
        Quantity traded =
            create_trade(inbound, current_order, MAX_QUANTITY, listener);
        if (traded > 0) {
          matched = true;
          // assert traded == current_quantity
          current_orders.erase(entry);
          inbound_qty -= traded;
        }
      } else {
        // current is AON, inbound is not AON.
        // inbound is not enough to satisfy current order's AON
        deferred_aons.push_back(entry);
      }
    } else {
      // neither are AON
      Quantity traded =
          create_trade(inbound, current_order, MAX_QUANTITY, listener);
      if (traded > 0) {
        matched = true;
        if (current_order.filled()) {
          current_orders.erase(entry);
        }
        inbound_qty -= traded;
      }
    }
  }
  return matched;
}

template <typename OrderPtr, template <typename...> typename Multimap,
          typename Listener>
bool OrderBook<OrderPtr, Multimap>::match_aon_order(
    Tracker &inbound, Price inbound_price, TrackerMap &current_orders,
    DeferredMatches &deferred_aons, Listener &&listener) {
  bool matched = false;
  Quantity inbound_qty = inbound.open_qty();
  Quantity deferred_qty = 0;

  DeferredMatches deferred_matches;

  typename TrackerMap::iterator pos = current_orders.begin();
  while (pos != current_orders.end() && !inbound.filled()) {
    auto entry = pos++;
    const ComparablePrice current_price = entry->first;
    if (!current_price.matches(inbound_price)) {
      // no more trades against current orders are possible
      break;
    }

    //////////////////////////////////////
    // Current price matches inbound price
    Tracker &current_order = entry->second;
    Quantity current_quantity = current_order.open_qty();

    if (current_order.all_or_none()) {
      // AON::AON
      // if the inbound order can satisfy the current order's AON condition
      if (current_quantity <= inbound_qty) {
        // if the the matched quantity can satisfy
        // the inbound order's AON condition
        if (inbound_qty <= current_quantity + deferred_qty) {
          // Try to create the deferred trades (if any) before creating
          // the trade with the current order.
          // What quantity will we need from the deferred orders?
          Quantity maxQty = inbound_qty - current_quantity;
          if (maxQty == try_create_deferred_trades(inbound, deferred_matches,
                                                   maxQty, maxQty,
                                                   current_orders, listener)) {
            inbound_qty -= maxQty;
            // finally execute this trade
            Quantity traded =
                create_trade(inbound, current_order, MAX_QUANTITY, listener);
            if (traded > 0) {
              // assert traded == current_quantity
              inbound_qty -= traded;
              matched = true;
              current_orders.erase(entry);
            }
          }
        } else {
          // AON::AON -- inbound could satisfy current, but
          // current cannot satisfy inbound;
          deferred_qty += current_quantity;
          deferred_matches.push_back(entry);
        }
      } else {
        // AON::AON -- inbound cannot satisfy current's AON
        deferred_aons.push_back(entry);
      }
    } else {
      // AON::REG

      // if we have enough to satisfy inbound
      if (inbound_qty <= current_quantity + deferred_qty) {
        Quantity traded = try_create_deferred_trades(
            inbound, deferred_matches,
            inbound_qty, // create as many as possible
            (inbound_qty > current_quantity)
                ? (inbound_qty - current_quantity)
                : 0, // but we need at least this many
            current_orders, listener);
        if (inbound_qty <= current_quantity + traded) {
          traded +=
              create_trade(inbound, current_order, MAX_QUANTITY, listener);
          if (traded > 0) {
            inbound_qty -= traded;
            matched = true;
          }
          if (current_order.filled()) {
            current_orders.erase(entry);
          }
        }
      } else {
        // not enough to satisfy inbound, yet.
        // remember the current order for later use
        deferred_qty += current_quantity;
        deferred_matches.push_back(entry);
      }
    }
  }
  return matched;
}
namespace {
const size_t AON_LIMIT = 5;
}

template <typename OrderPtr, template <typename...> typename Multimap,
          typename Listener>
Quantity OrderBook<OrderPtr, Multimap>::try_create_deferred_trades(
    Tracker &inbound, DeferredMatches &deferred_matches,
    Quantity maxQty, // do not exceed
    Quantity minQty, // must be at least
    TrackerMap &current_orders, Listener &&listener) {
  Quantity traded = 0;
  // create a vector of proposed trade quantities:
  std::vector<int> fills(deferred_matches.size());
  std::fill(fills.begin(), fills.end(), 0);
  Quantity foundQty = 0;
  auto pos = deferred_matches.begin();
  for (size_t index = 0; foundQty < maxQty && pos != deferred_matches.end();
       ++index) {
    auto entry = *pos++;
    Tracker &tracker = entry->second;
    Quantity qty = tracker.open_qty();
    // if this would put us over the limit
    if (foundQty + qty > maxQty) {
      if (tracker.all_or_none()) {
        qty = 0;
      } else {
        qty = maxQty - foundQty;
        // assert qty <= tracker.open_qty();
      }
    }
    foundQty += qty;
    fills[index] = qty;
  }

  if (foundQty >= minQty && foundQty <= maxQty) {
    // pass through deferred matches again, doing the trades.
    auto pos = deferred_matches.begin();
    for (size_t index = 0; traded < foundQty && pos != deferred_matches.end();
         ++index) {
      auto entry = *pos++;
      Tracker &tracker = entry->second;
      traded += create_trade(inbound, tracker, fills[index], listener);
      if (tracker.filled()) {
        current_orders.erase(entry);
      }
    }
  }
  return traded;
}

template <typename OrderPtr, template <typename...> typename Multimap,
          typename Listener>
Quantity OrderBook<OrderPtr, Multimap>::create_trade(Tracker &inbound_tracker,
                                                     Tracker &current_tracker,
                                                     Quantity maxQuantity,
                                                     Listener &&listener) {
  Price cross_price = current_tracker.ptr()->price();
  // If current order is a market order, cross at inbound price
  if (MARKET_ORDER_PRICE == cross_price) {
    cross_price = inbound_tracker.ptr()->price();
    if (MARKET_ORDER_PRICE == cross_price) {
      cross_price = marketPrice_;
      if (MARKET_ORDER_PRICE == cross_price) {
        // No price available for this order
        return 0;
      }
    }
  }
  Quantity fill_qty =
      std::min(maxQuantity, std::min(inbound_tracker.open_qty(),
                                     current_tracker.open_qty()));
  if (fill_qty > 0) {
    inbound_tracker.fill(fill_qty);
    current_tracker.fill(fill_qty);
    set_market_price(cross_price);

    listener(FillTag{}, inbound_tracker.ptr(), current_tracker.ptr(), fill_qty,
             cross_price, !inbound_tracker.open_qty(),
             !current_tracker.open_qty());
  }
  return fill_qty;
}

template <typename OrderPtr, template <typename...> typename Multimap>
void OrderBook<OrderPtr, Multimap>::move_callbacks(Callbacks &target) {
  COMPLAIN_ONCE("Ignoring call to deprecated method: move_callbacks");
  // We get to decide when callbacks happen.
  // And it *certainly* doesn't happen on another thread!
}

template <typename OrderPtr, template <typename...> typename Multimap>
void OrderBook<OrderPtr, Multimap>::perform_callbacks() {
  COMPLAIN_ONCE("Ignoring call to deprecated method: perform_callbacks");
  // We get to decide when callbacks happen.
}

// template <typename OrderPtr, template <typename...> typename Multimap>
// void OrderBook<OrderPtr, Multimap>::callback_now()
// {
//   for (auto &&cb : callbacks_)
//   {
//     perform_callback(cb);
//   }
//   callbacks_ = {};
// }

template <typename OrderPtr, template <typename...> typename Multimap>
std::ostream &OrderBook<OrderPtr, Multimap>::log(std::ostream &out) const {
  for (auto ask = asks_.rbegin(); ask != asks_.rend(); ++ask) {
    out << "  Ask " << ask->second.open_qty() << " @ " << ask->first
        << std::endl;
  }

  for (auto bid = bids_.begin(); bid != bids_.end(); ++bid) {
    out << "  Bid " << bid->second.open_qty() << " @ " << bid->first
        << std::endl;
  }
  return out;
}

} // namespace book
} // namespace liquibook
